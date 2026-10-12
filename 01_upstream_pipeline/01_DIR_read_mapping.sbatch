#!/bin/bash
#SBATCH -J DIR_pipeline
#SBATCH -c 24
#SBATCH --mem 128G
#SBATCH -t 72:00:00
#SBATCH -o slurm_logs/DIR_%j.out
#SBATCH -e slurm_logs/DIR_%j.err

# ============================================================
# 01_DIR_pipeline.sbatch  （v2）
# DIR 轨道（直接比对法）
# 流程：clean reads → BWA比对植物库 → idxstats/coverage/flagstat → 物种计数
#
# 【相对 v1 的改动】
#   1. 新增 flagstat.txt —— 给出比对率（v1 拿不到，因为 BAM 跑完即删）
#   2. 日志实时打印比对率，不用事后翻文件
#   3. 统计文件齐全后才删 BAM（防止中途失败丢数据）
#   4. 输出目录按批次隔离
# ============================================================
set -euo pipefail

trap 'echo -e "\n\033[31m===== DIR失败 =====\n样本：${SAMPLE:-未知}\n行号：$LINENO\n命令：$BASH_COMMAND\n退出码：$?\033[0m"' ERR

mkdir -p slurm_logs

##############################
## 1. 环境
##############################
echo "[$(date '+%H:%M:%S')] 加载环境..."
module load conda/2024.3
source $CONDA_ROOT/etc/profile.d/conda.sh
conda activate genome_index
echo "  ✓ 环境已激活"

##############################
## 2. 路径配置
##############################
CLEAN_DIR=${CLEAN_DIR:-~/accurate_clean_reads}

BATCH=${BATCH:-batch2}
OUT_BASE=${OUT_BASE:-~/plant_diet_results_${BATCH}}
mkdir -p $OUT_BASE

PLANT_BWA_IDX=~/seq_Table.fasta

THREADS=24
SORT_THREAD=$((THREADS / 2))

##############################
## 3. 索引校验
##############################
echo "[$(date '+%H:%M:%S')] 校验植物库BWA索引..."
for ext in amb ann bwt pac sa; do
    [ ! -f "${PLANT_BWA_IDX}.${ext}" ] && echo "  ✗ 缺失：${PLANT_BWA_IDX}.${ext}" && exit 1
done
echo "  ✓ BWA索引完整：$PLANT_BWA_IDX"

##############################
## 4. 批量处理
##############################
TOTAL=$(find $CLEAN_DIR -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) | wc -l)
echo "[$(date '+%H:%M:%S')] 共 $TOTAL 个样本，开始DIR轨道处理"

COUNT=0
for R1_FILE in $(find $CLEAN_DIR -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) | sort); do
    COUNT=$((COUNT + 1))

    BASENAME=$(basename "$R1_FILE")
    if [[ "$BASENAME" == *_R1.fq.gz ]]; then
        SAMPLE=$(basename "$R1_FILE" _R1.fq.gz)
        R2_FILE=${R1_FILE/_R1.fq.gz/_R2.fq.gz}
    else
        SAMPLE=$(basename "$R1_FILE" _1.fq.gz)
        R2_FILE=${R1_FILE/_1.fq.gz/_2.fq.gz}
    fi

    echo -e "\n[$(date '+%H:%M:%S')] [$COUNT/$TOTAL] $SAMPLE"
    [ ! -f "$R2_FILE" ] && echo "  ⚠ R2不存在，跳过" && continue

    SAMPLE_DIR=$OUT_BASE/$SAMPLE
    mkdir -p $SAMPLE_DIR

    # 断点续跑：★ v2 以 flagstat 为准（v1 的 idxstats 存在但没有比对率）
    if [ -f "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" ] && \
       [ -f "$SAMPLE_DIR/${SAMPLE}_idxstats.txt" ] && \
       [ -f "$SAMPLE_DIR/${SAMPLE}_coverage.txt" ]; then
        MAPPED=$(grep -m1 ' mapped (' "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" | sed 's/.*(\([0-9.]*\)%.*/\1/' || echo "?")
        echo "  ✓ 结果已存在（比对率 ${MAPPED}%），跳过"
        continue
    fi

    SORTED_BAM=$SAMPLE_DIR/${SAMPLE}_plant_sorted.bam

    # ===== 步骤1：BWA比对 =====
    echo "  [1/4] BWA-MEM比对植物库（24线程）..."
    bwa mem -t $THREADS -M $PLANT_BWA_IDX "$R1_FILE" "$R2_FILE" | \
        samtools sort -@ $SORT_THREAD -m 4G -o "$SORTED_BAM" -
    samtools index "$SORTED_BAM"
    echo "    ✓ 比对完成，BAM：$(du -h $SORTED_BAM | cut -f1)"

    # ===== 步骤2：idxstats =====
    echo "  [2/4] samtools idxstats..."
    samtools idxstats "$SORTED_BAM" > "$SAMPLE_DIR/${SAMPLE}_idxstats.txt"

    # ===== 步骤3：coverage =====
    echo "  [3/4] samtools coverage..."
    samtools coverage "$SORTED_BAM" > "$SAMPLE_DIR/${SAMPLE}_coverage.txt"

    # ===== 步骤4：flagstat（★ v2 新增：给出比对率）=====
    echo "  [4/5] samtools flagstat..."
    samtools flagstat "$SORTED_BAM" > "$SAMPLE_DIR/${SAMPLE}_flagstat.txt"
    TOTAL_RD=$(grep -m1 'in total' "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" | awk '{print $1}')
    MAPPED_RD=$(grep -m1 ' mapped (' "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" | awk '{print $1}')
    MAPPED_PCT=$(grep -m1 ' mapped (' "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" | sed 's/.*(\([0-9.]*\)%.*/\1/')
    echo "    总 reads：${TOTAL_RD}  比对：${MAPPED_RD}  比对率：${MAPPED_PCT}%"

    # ===== 步骤5：物种计数排序 =====
    echo "  [5/5] 物种reads计数排序..."
    awk '{print $1"\t"$3}' "$SAMPLE_DIR/${SAMPLE}_idxstats.txt" | \
        sort -k2 -nr > "$SAMPLE_DIR/${SAMPLE}_species_count.txt"

    # ★ 四个统计文件齐全后才删 BAM
    if [ -s "$SAMPLE_DIR/${SAMPLE}_idxstats.txt" ] && \
       [ -s "$SAMPLE_DIR/${SAMPLE}_coverage.txt" ] && \
       [ -s "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" ] && \
       [ -s "$SAMPLE_DIR/${SAMPLE}_species_count.txt" ]; then
        rm -f "$SORTED_BAM" "${SORTED_BAM}.bai"
        echo "  ✓ 完成（4 个统计文件已保存，BAM 已清理）"
    else
        echo "  ⚠ 统计文件不完整，保留 BAM 以便排查"
    fi
done

# ===== 汇总比对率 =====
SUMMARY=$OUT_BASE/mapping_rate_summary.tsv
{
  echo -e "Sample\tTotal_reads\tMapped_reads\tMapped_pct"
  for d in $OUT_BASE/*/; do
      s=$(basename "$d")
      f=$d/${s}_flagstat.txt
      [ -f "$f" ] || continue
      t=$(grep -m1 'in total' "$f" | awk '{print $1}')
      m=$(grep -m1 ' mapped (' "$f" | awk '{print $1}')
      p=$(grep -m1 ' mapped (' "$f" | sed 's/.*(\([0-9.]*\)%.*/\1/')
      echo -e "${s}\t${t}\t${m}\t${p}"
  done
} > "$SUMMARY"

echo -e "\n============================================================"
echo "[$(date '+%H:%M:%S')] DIR轨道全部完成！"
echo "  输出：$OUT_BASE"
echo "  每样本：idxstats.txt / coverage.txt / flagstat.txt / species_count.txt"
echo "  ★ 比对率汇总：$SUMMARY"
echo "============================================================"
awk -F'\t' 'NR>1{s+=$4; n++} END{if(n) printf "  平均比对率：%.2f%%（%d 个样本）\n", s/n, n}' "$SUMMARY"
