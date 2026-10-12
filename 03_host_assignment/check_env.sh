#!/usr/bin/env bash
# =============================================================================
#  00_check_env.sh —— 环境自检（不提交作业，只检查）
#
#  用法：  bash 00_check_env.sh
#
#  在正式提交 slurm 作业之前先跑这个，确认：
#    · conda 环境能激活
#    · bwa / samtools / python3 存在
#    · clean reads 找得到、R1/R2 配对完整
#    · 参考与位点文件完整，BWA 索引可建立
#  全部通过再去 sbatch。
# =============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CLEAN_DIR="${CLEAN_DIR:-$HOME/260412}"
REF_COMBINED="$SCRIPT_DIR/refs/muntjac_mito_combined.fna"

PASS=0; FAIL=0
ok(){ echo "  ✓ $*"; PASS=$((PASS+1)); }
no(){ echo "  ✗ $*"; FAIL=$((FAIL+1)); }
warn(){ echo "  ⚠ $*"; }

echo "============================================================"
echo " 环境自检"
echo "============================================================"
echo "脚本目录 : $SCRIPT_DIR"
echo "reads目录: $CLEAN_DIR"
echo

# ---------- 1. conda ----------
echo "[1/6] conda 环境"
if command -v conda >/dev/null 2>&1; then
    ok "conda 已在 PATH"
else
    if [ -n "${CONDA_ROOT:-}" ] && [ -f "$CONDA_ROOT/etc/profile.d/conda.sh" ]; then
        source "$CONDA_ROOT/etc/profile.d/conda.sh"
        ok "已从 CONDA_ROOT=$CONDA_ROOT 加载 conda"
    else
        no "找不到 conda，也不存在 \$CONDA_ROOT/etc/profile.d/conda.sh"
        warn "请在脚本里确认 module load conda/2024.3 那一行是否适用于你的集群"
    fi
fi

# ---------- 2. bwa / samtools / python3 ----------
echo "[2/6] 依赖程序"
for c in bwa samtools python3; do
    if command -v $c >/dev/null 2>&1; then
        ok "$c → $(command -v $c)"
    else
        no "$c 未找到（是否忘了 conda activate genome_index？）"
    fi
done

# ---------- 3. 附属文件 ----------
echo "[3/6] 附属文件"
for f in "$REF_COMBINED" "$SCRIPT_DIR/sites_all.tsv" "$SCRIPT_DIR/sites_all.bed" \
         "$SCRIPT_DIR/sites_transversion_only.tsv"; do
    if [ -f "$f" ]; then
        ok "$(basename "$f")  ($(stat -c%s "$f") 字节)"
    else
        no "缺失：$f"
    fi
done
if [ -f "$SCRIPT_DIR/sites_all.tsv" ]; then
    n=$(($(wc -l < "$SCRIPT_DIR/sites_all.tsv") - 1))
    if [ "$n" -eq 1894 ]; then
        ok "位点表行数正确：$n（947 位点 × 2 条染色体）"
    else
        no "位点表行数为 $n，应为 1894 —— 文件可能损坏"
    fi
fi

# ---------- 4. BWA 索引 ----------
echo "[4/6] BWA 索引（参考仅 16.4 kb，秒级）"
if command -v bwa >/dev/null 2>&1 && [ -f "$REF_COMBINED" ]; then
    if [ -f "${REF_COMBINED}.bwt" ]; then
        ok "索引已存在"
    else
        echo "     正在建立索引 ..."
        if bwa index "$REF_COMBINED" 2>/dev/null; then
            ok "索引建立成功"
        else
            no "索引建立失败"
        fi
    fi
    nref=$(grep -c '^>' "$REF_COMBINED")
    [ "$nref" -eq 2 ] && ok "参考含 $nref 条序列" || no "参考含 $nref 条序列，应为 2"
else
    warn "跳过（bwa 或参考不可用）"
fi

# ---------- 5. reads ----------
echo "[5/6] clean reads（$CLEAN_DIR）"
if [ ! -d "$CLEAN_DIR" ]; then
    no "目录不存在：$CLEAN_DIR"
    warn "用 pwd 确认当前家目录，用 ls 找到存放 fq.gz 的目录名，"
    warn "然后修改 05_mito_hostID.sbatch 里的 CLEAN_DIR= 这一行"
else
    n1=$(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*_1.fq.gz" | wc -l)
    nR1=$(find "$CLEAN_DIR" -maxdepth 1 -type f -name "*_R1.fq.gz" | wc -l)
    ok "找到 _1.fq.gz  $n1 个"
    [ "$nR1" -gt 0 ] && ok "找到 _R1.fq.gz $nR1 个"
    TOTAL=$((n1+nR1))
    [ "$TOTAL" -gt 0 ] && ok "样本总数：$TOTAL" || no "没有找到任何 R1 文件"

    # R1/R2 配对抽查
    miss=0; checked=0
    for R1 in $(find "$CLEAN_DIR" -maxdepth 1 -type f \( -name "*_1.fq.gz" -o -name "*_R1.fq.gz" \) | sort | head -5); do
        checked=$((checked+1))
        if [[ "$R1" == *_R1.fq.gz ]]; then R2=${R1/_R1.fq.gz/_R2.fq.gz}; else R2=${R1/_1.fq.gz/_2.fq.gz}; fi
        [ -f "$R2" ] || { no "缺配对：$(basename "$R1") → $(basename "$R2")"; miss=$((miss+1)); }
    done
    [ "$miss" -eq 0 ] && [ "$checked" -gt 0 ] && ok "抽查 $checked 个样本，R1/R2 配对正常"
fi

# ---------- 6. 磁盘 ----------
echo "[6/6] 磁盘空间"
AVAIL=$(df -BG "$HOME" | awk 'NR==2{gsub("G","",$4); print $4}')
if [ -n "$AVAIL" ]; then
    if [ "$AVAIL" -gt 20 ]; then
        ok "家目录可用 ${AVAIL} GB（输出约需 1-2 GB，充足）"
    else
        warn "家目录仅剩 ${AVAIL} GB，建议清理"
    fi
fi

echo
echo "============================================================"
echo " 结果：通过 $PASS 项，失败 $FAIL 项"
echo "============================================================"
if [ "$FAIL" -eq 0 ]; then
    echo " ✓ 一切正常，可以提交作业："
    echo "     cd $SCRIPT_DIR"
    echo "     mkdir -p slurm_logs"
    echo "     sbatch 05_mito_hostID.sbatch"
else
    echo " ✗ 请先解决上面标 ✗ 的问题，再提交作业。"
fi
echo
