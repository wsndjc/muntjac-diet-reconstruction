#!/usr/bin/env bash
#SBATCH -J mito_hostID
#SBATCH -c 24
#SBATCH --mem 32G
#SBATCH -t 24:00:00
#SBATCH -o slurm_logs/mito_%j.out
#SBATCH -e slurm_logs/mito_%j.err

# =============================================================================
# 05_mito_hostID.sbatch
# 线粒体宿主角逐判定 —— 「read 碱基投票」法
# 流程：clean reads → BWA比对两物种线粒体 → 判别位点碱基投票 → 二项检验
#
# 【为什么不用全基因组映射率】
#   原方法比较全基因组(2.49 Gb)映射率判定宿主，但两个参考组装差异极大：
#     染色体数     小麂 2n=46(23条)  vs  黑麂 2n=4(4条)   ← 麂属极端染色体重排
#     Scaffold N50  113,316,579       vs  646,574,120    ← 差 5.71 倍
#   实测：判定裕度中位数仅 0.02 个百分点，80% 样本裕度 <0.05 pp，
#         7 个样本两率完全相同，71% 样本倒向小麂一侧（正常应各半）。
#   文献：Atağ et al. 2023, MTaxi, Open Research Europe 2:100
#        "taxon identification using only the number of aligned bases in
#         comparative alignment is NOT RELIABLE"（绵羊/山羊 9 例判错 5 例，归因 N50）
#
# 【本方法原理】
#   不数 reads 总数，读每条 read 在种间差异位点上的碱基。
#   150 bp read 平均覆盖 150 × 5.79% ≈ 8.7 个差异位点；
#   来自小麂的 read 对小麂参考 ~0.15 错配(仅测序误差)，对黑麂参考 ~8.7 错配。
#   分离度约 58 倍 —— 即便发生交叉比对，read 自身碱基仍然可靠。
#
# 【资源】参考仅 16.4 kb，比对开销可忽略；32G 内存绰绰有余
# =============================================================================

set -euo pipefail

trap 'echo -e "\n\033[31m===== 失败 =====\n样本：${SAMPLE:-未知}\n行号：$LINENO\n命令：$BASH_COMMAND\n\033[0m"' ERR

mkdir -p slurm_logs

##############################
## 1. 环境
##############################
echo "[$(date '+%H:%M:%S')] 加载环境..."
if command -v conda >/dev/null 2>&1; then
    echo "  ✓ conda 已在 PATH"
else
    module load conda/2024.3 2>/dev/null || true
    if [ -n "${CONDA_ROOT:-}" ] && [ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]; then
        source "$CONDA_ROOT/etc/profile.d/conda.sh"
    fi
fi
conda activate genome_index || { echo "✗ 无法激活 genome_index 环境"; exit 1; }
echo "  ✓ 环境已激活：$CONDA_DEFAULT_ENV"
for c in bwa samtools python3; do
    command -v $c >/dev/null 2>&1 || { echo "✗ 缺少 $c"; exit 1; }
done
echo "  ✓ bwa / samtools / python3 就绪"

##############################
## 2. 路径配置
##############################
CLEAN_DIR=${CLEAN_DIR:-$HOME/260412}          # 双端 clean reads
OUT_BASE=${OUT_BASE:-$HOME/mito_hostID_results}
# 自动识别本脚本所在目录（sbatch 从哪个目录提交就指哪个）
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(pwd)}"

if [ ! -f "$SCRIPT_DIR/refs/muntjac_mito_combined.fna" ]; then
    echo "✗ 在 $SCRIPT_DIR 找不到 refs/muntjac_mito_combined.fna"
    echo "  请先 cd 到解压后的 mito_id_bundle 目录，再执行 sbatch"
    exit 1
fi
echo "  脚本目录：$SCRIPT_DIR"
echo "  reads目录：$CLEAN_DIR"

REF_COMBINED=$SCRIPT_DIR/refs/muntjac_mito_combined.fna
SITES_ALL=$SCRIPT_DIR/sites_all.tsv
SITES_TV=$SCRIPT_DIR/sites_transversion_only.tsv

THREADS=24
# 并行策略：线粒体参考仅 16.4 kb，比对几乎不耗时，瓶颈是 bwa 内部【单线程】gzip 解压。
# 若一次只跑一个样本并给满 24 线程，则 23 个核全程闲置在等解压。
# 改为 N 个样本并行、每样本分到 THREADS/N 线程，可近线性加速。
PARALLEL=${PARALLEL:-6}
PER_THREADS=$(( THREADS / PARALLEL ))
[ "$PER_THREADS" -lt 1 ] && PER_THREADS=1
SORT_THREAD=2
MIN_MAPQ=20
MIN_BQ=20

mkdir -p $OUT_BASE/{bam,log,votes}

##############################
## 3. 校验
##############################
echo "[$(date '+%H:%M:%S')] 校验输入与索引..."
for f in "$REF_COMBINED" "$SITES_ALL" "$SITES_TV"; do
    [ -f "$f" ] || { echo "  ✗ 缺失：$f"; exit 1; }
done
echo "  ✓ 附属文件完整"

if [ ! -f "${REF_COMBINED}.bwt" ]; then
    echo "  建立BWA索引（16 kb，秒级）..."
    bwa index "$REF_COMBINED"
fi
for ext in amb ann bwt pac sa; do
    [ ! -f "${REF_COMBINED}.${ext}" ] && echo "  ✗ 索引缺失：${ext}" && exit 1
done
echo "  ✓ BWA索引完整"

##############################
## 4. 投票脚本
##############################
cat > $OUT_BASE/vote_sites.py <<'PYEOF'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""从 samtools mpileup(-s) 统计判别位点碱基，做二项检验判定宿主物种。

以「read 自身携带的碱基」投票，而非数 reads 总数。
即使 read 同时能比上两个参考（交叉比对），其碱基仍然指向真实来源物种。

注意：合并参考含两条染色体，同一分化事件在两个坐标系各出现一次
（chr A 的 posA 与 chr B 的 posB）。位点表统一给出
(base_reevesi, base_crinifrons) 两列，因此无论 read 落在哪条染色体上，
其碱基都能被正确归类到物种，不会互相抵消。
"""
import sys
from math import lgamma, log, exp
from collections import Counter

def binom_two_sided(k, n, p=0.5):
    """精确二项检验（双侧），对数域避免溢出。"""
    if n <= 0:
        return 1.0
    def lpmf(i):
        return (lgamma(n + 1) - lgamma(i + 1) - lgamma(n - i + 1)
                + i * log(p) + (n - i) * log(1 - p))
    obs = lpmf(k)
    tot = 0.0
    for i in range(n + 1):
        v = lpmf(i)
        if v <= obs + 1e-9:
            tot += exp(v)
    return min(1.0, tot)

def main():
    sites_file, mpileup_file, min_bq, sample = sys.argv[1:5]
    min_bq = int(min_bq)

    # chrom:pos -> (base_reevesi, base_crinifrons)
    sites = {}
    with open(sites_file, encoding="utf-8") as f:
        for line in f:
            if line.startswith("chrom") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                continue
            sites[(c[0], int(c[1]))] = (c[2].upper(), c[3].upper())

    n_ree = n_cri = n_amb = 0
    per_site = Counter()
    with open(mpileup_file, encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 5:
                continue
            key = (c[0], int(c[1]))
            if key not in sites:
                continue
            base_ree, base_cri = sites[key]
            ref = c[2].upper()
            rb = c[4]
            # 逐碱基质量（-s 输出的第 6 列，空格分隔）
            quals = []
            if len(c) > 5:
                for tok in c[5].split():
                    try:
                        quals.append(int(tok))
                    except ValueError:
                        pass
            qi = 0
            i = 0
            while i < len(rb):
                ch = rb[i]
                if ch == "^":            # read 起始标记，后跟一个 mapq 字符
                    i += 2
                    continue
                if ch == "$":
                    i += 1
                    continue
                if ch in "+-":           # indel 指示：+2AG / -3TTC
                    i += 1
                    num = ""
                    while i < len(rb) and rb[i].isdigit():
                        num += rb[i]; i += 1
                    i += int(num) if num else 0
                    continue
                if ch in "*#":           # 缺失
                    i += 1; qi += 1
                    continue
                if ch in ".,":
                    b = ref
                elif ch.upper() in "ACGTN":
                    b = ch.upper()
                else:
                    i += 1; qi += 1
                    continue
                bq = quals[qi] if qi < len(quals) else 30
                qi += 1
                i += 1
                if bq < min_bq or b == "N":
                    continue
                if b == base_ree:
                    n_ree += 1; per_site[key] += 1
                elif b == base_cri:
                    n_cri += 1; per_site[key] += 1
                else:
                    n_amb += 1

    n = n_ree + n_cri
    p = binom_two_sided(n_cri, n) if n else 1.0
    if n == 0:
        call = "NA"
    elif p < 0.05 and n_ree > n_cri:
        call = "reevesi"
    elif p < 0.05 and n_cri > n_ree:
        call = "crinifrons"
    else:
        call = "ambiguous"
    print(f"{sample}\t{n_ree}\t{n_cri}\t{n_amb}\t{n}\t{len(per_site)}\t{p:.6g}\t{call}")

if __name__ == "__main__":
    main()
PYEOF
chmod +x $OUT_BASE/vote_sites.py

##############################
## 5. 批量处理
##############################
TOTAL=$(find $CLEAN_DIR -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) | wc -l)
echo "[$(date '+%H:%M:%S')] 共 $TOTAL 个样本，开始线粒体宿主判定"

# 结果表头
RESULT_ALL=$OUT_BASE/host_call_allsites.tsv
RESULT_TV=$OUT_BASE/host_call_transversion.tsv
echo -e "SampleID\tvotes_reevesi\tvotes_crinifrons\tambiguous\tn_votes\tsites_hit\tp_value\tcall" > $RESULT_ALL
echo -e "SampleID\tvotes_reevesi\tvotes_crinifrons\tambiguous\tn_votes\tsites_hit\tp_value\tcall" > $RESULT_TV

##############################
## 5. 单样本处理函数
##############################
run_one () {
    local R1_FILE="$1" IDX="$2" TOTAL="$3"
    local BASENAME SAMPLE R2_FILE BAM CALL_A CALL_T
    BASENAME=$(basename "$R1_FILE")
    if [[ "$BASENAME" == *_R1.fq.gz ]]; then
        SAMPLE=$(basename "$R1_FILE" _R1.fq.gz)
        R2_FILE=${R1_FILE/_R1.fq.gz/_R2.fq.gz}
    else
        SAMPLE=$(basename "$R1_FILE" _1.fq.gz)
        R2_FILE=${R1_FILE/_1.fq.gz/_2.fq.gz}
    fi

    if [ ! -f "$R2_FILE" ]; then
        echo "[$IDX/$TOTAL] ⚠ $SAMPLE R2不存在，跳过"; return
    fi

    # 断点续跑
    if [ -f "$OUT_BASE/log/${SAMPLE}.mpileup" ]; then
        echo "[$IDX/$TOTAL] ✓ $SAMPLE mpileup已存在，跳过比对"
    else
        BAM=$OUT_BASE/bam/${SAMPLE}_mito.bam
        echo "[$(date '+%H:%M:%S')] [$IDX/$TOTAL] $SAMPLE 比对开始（${PER_THREADS}线程）"
        # 关键优化：-F 4 在 sort 之前丢弃未比对 reads。
        # bwa mem 默认把未比对 reads 也写进 SAM，114M reads 会产生 ~40 GB SAM，
        # 而能落在 947 个判别位点上的只有几百到几千条。
        if bwa mem -t $PER_THREADS -M "$REF_COMBINED" "$R1_FILE" "$R2_FILE" \
                2>$OUT_BASE/log/${SAMPLE}.bwa.log \
             | samtools view -@ 1 -b -F 4 - \
             | samtools sort -@ $SORT_THREAD -m 1G -o "$BAM" - ; then
            samtools index "$BAM"
            echo "    ✓ $SAMPLE 命中线粒体 $(samtools view -c $BAM) reads"
        else
            echo "    ✗ $SAMPLE 比对失败，见 log/${SAMPLE}.bwa.log"; return
        fi

        samtools mpileup -f "$REF_COMBINED" -l $SCRIPT_DIR/sites_all.bed \
            -q $MIN_MAPQ -Q $MIN_BQ -s -d 100000 -B "$BAM" 2>/dev/null \
            > $OUT_BASE/log/${SAMPLE}.mpileup
        rm -f "$BAM" "${BAM}.bai"
    fi

    # 投票：写各自独立文件，避免并行追加冲突（最后再合并）
    python3 $OUT_BASE/vote_sites.py "$SITES_ALL" $OUT_BASE/log/${SAMPLE}.mpileup $MIN_BQ "$SAMPLE" \
        > $OUT_BASE/votes/${SAMPLE}.all.tsv
    python3 $OUT_BASE/vote_sites.py "$SITES_TV"  $OUT_BASE/log/${SAMPLE}.mpileup $MIN_BQ "$SAMPLE" \
        > $OUT_BASE/votes/${SAMPLE}.tv.tsv
    CALL_A=$(cut -f8 $OUT_BASE/votes/${SAMPLE}.all.tsv)
    CALL_T=$(cut -f8 $OUT_BASE/votes/${SAMPLE}.tv.tsv)
    echo "[$(date '+%H:%M:%S')] [$IDX/$TOTAL] ✓ $SAMPLE 完成｜全部位点→${CALL_A}｜仅颠换→${CALL_T}"
}

export -f run_one
export OUT_BASE REF_COMBINED SCRIPT_DIR SITES_ALL SITES_TV
export MIN_MAPQ MIN_BQ PER_THREADS SORT_THREAD

##############################
## 6. 并行调度
##############################
find $CLEAN_DIR -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) | sort | \
    awk -v tot="$TOTAL" '{printf "%s|%d|%d\n", $0, NR, tot}' | \
    xargs -P "$PARALLEL" -I{} bash -c 'IFS="|" read -r f i t <<< "{}"; run_one "$f" "$i" "$t"'

##############################
## 7. 合并结果
##############################
for tag in all tv; do
    if [ "$tag" = "all" ]; then OUT_F=$RESULT_ALL; else OUT_F=$RESULT_TV; fi
    { echo -e "SampleID\tvotes_reevesi\tvotes_crinifrons\tambiguous\tn_votes\tsites_hit\tp_value\tcall"
      cat $OUT_BASE/votes/*.$tag.tsv 2>/dev/null | sort -k1,1V; } > $OUT_F
done

echo -e "\n============================================================"
echo "[$(date '+%H:%M:%S')] 全部完成！"
echo "  全部位点(947)：$RESULT_ALL"
echo "  仅颠换(83)   ：$RESULT_TV"
echo "============================================================"

python3 - "$RESULT_ALL" "$RESULT_TV" <<'PYEOF'
import sys
from collections import Counter
for path, label in [(sys.argv[1], "全部位点(947)"), (sys.argv[2], "仅颠换(83)")]:
    rows = []
    with open(path, encoding="utf-8") as f:
        next(f, None)
        for ln in f:
            c = ln.rstrip("\n").split("\t")
            if len(c) >= 8:
                rows.append(c)
    cnt = Counter(r[7] for r in rows)
    print(f"\n【{label}】共 {len(rows)} 个样本")
    print(f"  reevesi(小麂)     {cnt.get('reevesi', 0):>3}")
    print(f"  crinifrons(黑麂)  {cnt.get('crinifrons', 0):>3}")
    print(f"  ambiguous         {cnt.get('ambiguous', 0):>3}")
    print(f"  NA(无判别reads)   {cnt.get('NA', 0):>3}")
PYEOF
