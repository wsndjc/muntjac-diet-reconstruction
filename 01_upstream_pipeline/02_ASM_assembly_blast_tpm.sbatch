#!/bin/bash
#SBATCH -J ASM_batch2
#SBATCH -c 24
#SBATCH --mem 200G
#SBATCH -t 72:00:00
#SBATCH --array=1-200%3
#SBATCH -o slurm_logs/ASM_%A_%a.out
#SBATCH -e slurm_logs/ASM_%A_%a.err

# =============================================================================
# 02_ASM_pipeline.sbatch   （第二批 · 与第一批参数逐字一致）
#
# 流程：clean reads → MEGAHIT组装 → BLAST注释(植物库) → Bowtie2回比 → TPM
#
# 【与第一批保持一致的部分】—— 全部照搬 final_assembly_blast_csv.sh，未改一字
#   MIN_CONTIG_LEN=300
#   BLAST_EVALUE=1e-8
#   blastn -task blastn -outfmt "6 qseqid sseqid pident length ... sscinames"
#   blastn -max_target_seqs 1 -num_threads 24
#   megahit -m 0.5
#   bowtie2-build --threads 24
#   bowtie2 -p 24 -x <idx> -1 <R1> -2 <R2> | samtools sort -@ 12
#   TPM 的 awk 计算式（rpk = reads/(len/1000)，再归一化到 1e6）
#   ★ 不加 -F 4（第一批没有，保持一致）
#
# 【只改了资源调度，不影响科学结果】
#   1. 文件名匹配：第一批是 *_R1.fq.gz，第二批是 *_1.fq.gz
#      改为同时匹配两种（否则循环匹配 0 个文件，脚本空转）
#   2. 数组并行：串行 86 样本 × 2-4h = 172-344h 会超时
#      改为每样本一个作业，%3 同时跑 3 个
#   3. 内存：128G → 200G
#   4. 输出目录按批次隔离，不覆盖第一批
# =============================================================================

set -euo pipefail

trap 'echo -e "\n\033[31m===== ASM失败 =====\n样本：${SAMPLE:-未知}\n行号：$LINENO\n命令：$BASH_COMMAND\033[0m"' ERR

mkdir -p slurm_logs

##############################
## 1. 环境
##############################
echo "[$(date '+%H:%M:%S')] 加载环境..."
module load conda/2024.3
source $CONDA_ROOT/etc/profile.d/conda.sh
conda activate genome_index
echo "  ✓ 环境已激活"
which megahit blastn bowtie2 samtools

##############################
## 2. 路径配置
##############################
CLEAN_DIR=${CLEAN_DIR:-~/accurate_clean_reads}

# 输出目录按批次隔离
BATCH=${BATCH:-batch2}
OUT_BASE=${OUT_BASE:-~/assembly_blast_mapping_${BATCH}}

# 植物参考序列与BLAST数据库（与第一批同一个库）
PLANT_FA=~/seq_Table.fasta
BLAST_DB_DIR=~/plant_blastdb
BLAST_DB=$BLAST_DB_DIR/seq_Table
mkdir -p $OUT_BASE $BLAST_DB_DIR

# 参数（与第一批逐字一致）
MIN_CONTIG_LEN=300
BLAST_EVALUE=1e-8
THREADS=24
MEGAHIT_MEM=0.5

##############################
## 3. 植物库BLAST索引（幂等，与第一批同一个库）
##############################
if [ ! -f "${BLAST_DB}.nhr" ]; then
    echo "[$(date '+%H:%M:%S')] 首次运行：为 seq_Table.fasta 构建BLAST数据库..."
    bwa_index_present=0
    [ ! -f "$PLANT_FA" ] && { echo "  ✗ 缺少 $PLANT_FA"; exit 1; }
    makeblastdb -dbtype nucl \
                -in $PLANT_FA \
                -input_type fasta \
                -parse_seqids \
                -out $BLAST_DB
    echo "  ✓ BLAST库构建完成"
else
    echo "  ✓ BLAST植物库已存在：$BLAST_DB"
fi

##############################
## 4. 取本任务负责的样本
##############################
SAMPLE_LIST=$OUT_BASE/sample_list.txt
if [ ! -f "$SAMPLE_LIST" ]; then
    find "$CLEAN_DIR" -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) \
        | sed 's|.*/||; s|_1\.fq\.gz$||; s|_R1\.fq\.gz$||' | sort -u > "$SAMPLE_LIST"
fi

N=$(wc -l < "$SAMPLE_LIST")
IDX=${SLURM_ARRAY_TASK_ID:-1}

echo "[$(date '+%H:%M:%S')] 任务 $IDX / 共 $N 个样本"
if [ "$IDX" -gt "$N" ]; then
    echo "  任务编号超出样本数，正常退出"
    exit 0
fi

SAMPLE=$(sed -n "${IDX}p" "$SAMPLE_LIST")
[ -z "$SAMPLE" ] && { echo "  ✗ 取不到样本名"; exit 1; }

# 自动兼容 _1.fq.gz / _R1.fq.gz 两种文件名
if   [ -f "${CLEAN_DIR}/${SAMPLE}_1.fq.gz" ];  then
    R1_FILE=${CLEAN_DIR}/${SAMPLE}_1.fq.gz;  R2_FILE=${CLEAN_DIR}/${SAMPLE}_2.fq.gz
elif [ -f "${CLEAN_DIR}/${SAMPLE}_R1.fq.gz" ]; then
    R1_FILE=${CLEAN_DIR}/${SAMPLE}_R1.fq.gz; R2_FILE=${CLEAN_DIR}/${SAMPLE}_R2.fq.gz
else
    echo "  ✗ 找不到 ${SAMPLE} 的 reads，跳过"; exit 0
fi

[ ! -f "$R2_FILE" ] && echo "  ⚠️ 跳过 $SAMPLE：R2文件不存在" && exit 0

SAMPLE_DIR=$OUT_BASE/$SAMPLE
mkdir -p $SAMPLE_DIR
echo -e "\n====================================="
echo "🎯 处理样本：$SAMPLE"
echo "   reads：$(basename $R1_FILE) + $(basename $R2_FILE)"
echo "====================================="

##############################
## 5. 步骤1：MEGAHIT 组装 Contigs
##############################
CONTIG_FA=$SAMPLE_DIR/${SAMPLE}_contigs.fa
if [ ! -s "$CONTIG_FA" ]; then
    echo "📦 1. 组装中：megahit 最小长度${MIN_CONTIG_LEN}"
    megahit --min-contig-len $MIN_CONTIG_LEN \
            -t $THREADS \
            -m $MEGAHIT_MEM \
            -1 "$R1_FILE" \
            -2 "$R2_FILE" \
            -o "$SAMPLE_DIR/megahit_tmp" \
            > "$SAMPLE_DIR/${SAMPLE}.megahit.log" 2>&1

    cp "$SAMPLE_DIR/megahit_tmp/final.contigs.fa" "$CONTIG_FA"
    rm -rf "$SAMPLE_DIR/megahit_tmp"
    echo "  ✅ 组装完成：$(grep -c '^>' $CONTIG_FA) 条 contig"
else
    echo "  ✅ 已存在组装结果，跳过组装（$(grep -c '^>' $CONTIG_FA) 条 contig）"
fi

# 空contig检测
[ ! -s "$CONTIG_FA" ] && echo "  ❌ $SAMPLE 无有效contig，跳过后续步骤" && exit 0

##############################
## 6. 步骤2：BLAST 物种注释（植物库）
##############################
BLAST_OUT=$SAMPLE_DIR/${SAMPLE}_blast_annotation.csv
if [ ! -s "$BLAST_OUT" ]; then
    echo "🔍 2. BLAST 物种注释"
    blastn -query "$CONTIG_FA" \
           -db "$BLAST_DB" \
           -task blastn \
           -evalue $BLAST_EVALUE \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames" \
           -max_target_seqs 1 \
           -num_threads $THREADS | tr '\t' ',' > "$BLAST_OUT.body"

    # 为 CSV 文件添加表头
    echo "qseqid,sseqid,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,staxids,sscinames" \
        | cat - "$BLAST_OUT.body" > "$BLAST_OUT"
    rm -f "$BLAST_OUT.body"

    N_HIT=$(($(wc -l < "$BLAST_OUT") - 1))
    echo "  ✅ 注释完成：${N_HIT} 条 contig 有植物比对"
    if [ "$N_HIT" -eq 0 ]; then
        echo "  ⚠⚠ 没有任何 contig 比对上植物库！请检查 $BLAST_DB"
    fi
else
    echo "  ✅ 注释结果已存在，跳过"
fi

##############################
## 7. 步骤3：Bowtie2 回比 + TPM 丰度计算
##############################
echo "📊 3. Bowtie2 回比并计算 TPM"
BAM_FILE=$SAMPLE_DIR/${SAMPLE}_contig_mapping.bam
TPM_OUT=$SAMPLE_DIR/${SAMPLE}_tpm_abundance.csv

if [ ! -s "$TPM_OUT" ]; then
    # Bowtie2 建库
    bowtie2-build --threads $THREADS "$CONTIG_FA" "$SAMPLE_DIR/${SAMPLE}_bt2_index" \
        > "$SAMPLE_DIR/${SAMPLE}.bt2build.log" 2>&1

    # Bowtie2 比对 + 排序（与第一批一致，不加 -F 4）
    bowtie2 -p $THREADS -x "$SAMPLE_DIR/${SAMPLE}_bt2_index" \
            -1 "$R1_FILE" -2 "$R2_FILE" 2>"$SAMPLE_DIR/${SAMPLE}.bowtie2.log" | \
    samtools sort -@ 12 -o "$BAM_FILE" -
    samtools index "$BAM_FILE"

    # 计算 TPM 丰度，直接输出标准 CSV 格式（awk 逻辑与第一批逐字一致）
    samtools idxstats "$BAM_FILE" | awk -F '\t' '
    {
        if ($3 > 0) {
            rpk = $3 / ($2 / 1000);
            rpk_sum += rpk;
            rpk_dict[$1] = rpk;
            len_dict[$1] = $2;
            reads_dict[$1] = $3;
        }
    }
    END {
        print "Contig,Length,Mapped_Reads,TPM"
        for (c in rpk_dict) {
            tpm = (rpk_dict[c] / rpk_sum) * 1000000;
            printf "%s,%d,%d,%.2f\n", c, len_dict[c], reads_dict[c], tpm;
        }
    }' > "$TPM_OUT"

    # ★ 附加输出：回比率（第一批没有这一步，但保留 BAM 后可随时补算）
    samtools flagstat "$BAM_FILE" > "$SAMPLE_DIR/${SAMPLE}_flagstat.txt"
    MAPPED_PCT=$(grep -m1 ' mapped (' "$SAMPLE_DIR/${SAMPLE}_flagstat.txt" | sed 's/.*(\([0-9.]*\)%.*/\1/' || echo "?")
    echo "  ✅ TPM完成，定量contig：$(($(wc -l < "$TPM_OUT") - 1)) 条；回比率 ${MAPPED_PCT}%"

    # 清理临时索引（★ 保留 BAM，与第一批一致）
    rm -f "$SAMPLE_DIR"/${SAMPLE}_bt2_index*.bt2
else
    echo "  ✅ TPM结果已存在，跳过"
fi

echo -e "\n🎉 $SAMPLE 全流程完成！"
echo "  - 组装contig：$CONTIG_FA"
echo "  - BLAST注释：$BLAST_OUT"
echo "  - TPM丰度  ：$TPM_OUT"
