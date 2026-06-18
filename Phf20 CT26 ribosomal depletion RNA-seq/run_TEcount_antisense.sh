#!/usr/bin/env bash

# ==========================================
# TEcount antisense quantification
# Project: Phf20 CT26 ribosomal depletion RNA-seq
#
# The library is reverse-stranded / first-strand.
# Existing sense TEcount used: --stranded reverse
# Antisense TEcount should use: --stranded forward
# ==========================================

# 1. 初始化 conda
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    # shellcheck disable=SC1091
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

# 2. 激活环境
conda activate te_env

ulimit -n 65535

# [路径配置]
WORKDIR="/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts_antisense_TEcount"

# 小鼠注释
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"

# 并发数。TEcount 会占用内存，4 个样本并行通常足够。
MAX_JOBS=4

mkdir -p "${COUNTS_DIR}"

valid_bam() {
    local bam_file="$1"
    [ -s "$bam_file" ] && samtools quickcheck "$bam_file" >/dev/null 2>&1
}

echo "=== TEcount 反义链定量开始 ==="
echo "WORKDIR: ${WORKDIR}"
echo "ALIGN_DIR: ${ALIGN_DIR}"
echo "COUNTS_DIR: ${COUNTS_DIR}"
echo "GTF_GENE: ${GTF_GENE}"
echo "GTF_TE: ${GTF_TE}"
echo "Stranded mode: forward (antisense for reverse-stranded library)"

if ! command -v TEcount >/dev/null 2>&1; then
    echo "❌ [报错] 未找到 TEcount，请检查 te_env 环境。"
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "❌ [报错] 未找到 samtools，请检查 te_env 环境。"
    exit 1
fi

if [ ! -f "${GTF_GENE}" ]; then
    echo "❌ [报错] 基因 GTF 不存在: ${GTF_GENE}"
    exit 1
fi

if [ ! -f "${GTF_TE}" ]; then
    echo "❌ [报错] TE GTF 不存在: ${GTF_TE}"
    exit 1
fi

shopt -s nullglob
bam_files=("${ALIGN_DIR}"/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

if [ ${#bam_files[@]} -eq 0 ]; then
    echo "❌ [报错] 未找到 BAM 文件: ${ALIGN_DIR}/*.Aligned.sortedByCoord.out.bam"
    exit 1
fi

for bam_file in "${bam_files[@]}"; do
    sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
    project_prefix="${COUNTS_DIR}/${sample_name}_antisense"
    output_file="${project_prefix}.cntTable"
    log_file="${project_prefix}.TEcount.log"

    if ! valid_bam "${bam_file}"; then
        echo "⚠️ [跳过] ${sample_name} 不是有效 BAM: ${bam_file}"
        continue
    fi

    if [ -s "${output_file}" ]; then
        echo "✅ [跳过] ${sample_name} 反义链 TEcount 已完成: ${output_file}"
        continue
    fi

    echo ">>> 启动 ${sample_name} 反义链 TEcount"
    (
        TEcount --sortByPos --format BAM --mode multi \
                --GTF "${GTF_GENE}" --TE "${GTF_TE}" \
                --project "${project_prefix}" \
                --stranded forward -b "${bam_file}" \
                > "${log_file}" 2>&1 \
        && echo "🎉 [TEcount antisense 完成] ${sample_name}" \
        || echo "❌ [TEcount antisense 失败] ${sample_name}，详见 ${log_file}"
    ) &

    while (( $(jobs -r -p | wc -l) >= MAX_JOBS )); do
        sleep 10
    done
done

wait

echo "✅ === TEcount 反义链定量结束 ==="
echo "输出目录: ${COUNTS_DIR}"
