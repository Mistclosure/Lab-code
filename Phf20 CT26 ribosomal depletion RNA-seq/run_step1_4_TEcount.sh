#!/usr/bin/env bash

# ==========================================
# TE Analysis Pipeline Step 1-4
# Project: Phf20 CT26 ribosomal depletion RNA-seq
# Includes: fastp + bowtie2 rRNA removal + STAR alignment + TEcount
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

# ======================
# 核心配置区域
# ======================

# 解除 Linux 文件打开数量限制
ulimit -n 65535

# [CPU 策略]
HIGH_THREADS=100
MID_THREADS=80
LOW_THREADS=50

# [内存 策略]
DUMP_MEM="8000MB"
STAR_RAM="150000000000"

# [路径配置]
WORKDIR="/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
RAW_DIR="${WORKDIR}/rawdata"
TRIM_DIR="${WORKDIR}/trimmed_fastq"
CLEAN_DIR="${WORKDIR}/clean_non_rRNA"
ALIGN_DIR="${WORKDIR}/alignments"
COUNTS_DIR="${WORKDIR}/counts"
BOWTIE_TMP_DIR="/mnt/disk1/qiuzerui/expriments/Phf20_CT26_ribo_bowtie_tmp"

# 小鼠注释与索引
ANNO_DIR="/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38"
STAR_INDEX="/mnt/windowsdata/qiuzerui/indexes/star_index_m39"
GTF_GENE="${ANNO_DIR}/gencode.vM38.annotation_PRI.gtf"
GTF_TE="${ANNO_DIR}/m39_TE.gtf"
RRNA_INDEX="${ANNO_DIR}/rRNA_mtDNA_index"

# 初始化目录
echo ">>> 正在初始化目录..."
mkdir -p "${TRIM_DIR}" "${CLEAN_DIR}" "${ALIGN_DIR}" "${COUNTS_DIR}" "${BOWTIE_TMP_DIR}"

valid_bam() {
    local bam_file="$1"
    [ -s "$bam_file" ] && samtools quickcheck "$bam_file" >/dev/null 2>&1
}


# ==========================================
# Step 1-3: 高兼容性匹配模式
# 支持 .fastq / .fq / .fastq.gz / .fq.gz
# 当前项目 rawdata 文件示例:
# L1MLF0204314-Scr_1_Mixt.R1.raw.fastq
# ==========================================
echo "=== Step 1-3: 正在扫描 ${RAW_DIR} 下的 FASTQ 文件 ==="

shopt -s nullglob
all_files=(
    "${RAW_DIR}"/*.fastq
    "${RAW_DIR}"/*.fq
    "${RAW_DIR}"/*.fastq.gz
    "${RAW_DIR}"/*.fq.gz
)
shopt -u nullglob

if [ ${#all_files[@]} -gt 0 ]; then
    for r1_file in "${all_files[@]}"; do
        filename=$(basename "${r1_file}")

        # 只处理明确的 R1 文件，避免把样本名中的 _1_/_2_ 当成读段编号。
        if [[ "$filename" =~ ^(.+)([_.]R1)(.*\.(fastq|fq)(\.gz)?)$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            r1_id="${BASH_REMATCH[2]}"
            suffix="${BASH_REMATCH[3]}"
            r2_id="${r1_id/R1/R2}"
            r2_filename="${sample_name}${r2_id}${suffix}"
            r2_file="${RAW_DIR}/${r2_filename}"
        elif [[ "$filename" =~ ^(.+)([_.]1)(\.(fastq|fq)(\.gz)?)$ ]]; then
            sample_name="${BASH_REMATCH[1]}"
            r1_id="${BASH_REMATCH[2]}"
            suffix="${BASH_REMATCH[3]}"
            r2_id="${r1_id/1/2}"
            r2_filename="${sample_name}${r2_id}${suffix}"
            r2_file="${RAW_DIR}/${r2_filename}"
        else
            continue
        fi

        # 检查推导出的 R2 文件是否存在
        if [ ! -f "$r2_file" ]; then
            echo "⚠️ [警告] 发现 R1 但未找到对应 R2: $filename (预期 R2: $r2_filename)"
            continue
        fi

        # 检查断点续跑：比对是否已完成
        bam_out="${ALIGN_DIR}/${sample_name}.Aligned.sortedByCoord.out.bam"
        if valid_bam "${bam_out}"; then
            echo "✅ [跳过] ${sample_name} 已完成比对。"
            continue
        elif [ -f "${bam_out}" ]; then
            echo "⚠️ [重跑] ${sample_name} 已有 BAM 但无效或不完整。"
            rm -f "${bam_out}" "${bam_out}.bai"
        fi

        echo ">>> 🚀 正在处理样本: ${sample_name} <<<"

        # [1/3] Fastp 质控
        if [ ! -f "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ]; then
            fastp -i "${r1_file}" -I "${r2_file}" \
                  -o "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                  -O "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                  -h "${TRIM_DIR}/${sample_name}_fastp.html" \
                  -j "${TRIM_DIR}/${sample_name}_fastp.json" \
                  --thread ${LOW_THREADS} --detect_adapter_for_pe --length_required 25 2> /dev/null
        fi

        if [ ! -s "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" ] || [ ! -s "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" ]; then
            echo "❌ [报错] ${sample_name} fastp 输出缺失，跳过该样本。"
            continue
        fi

        # [2/3] Bowtie2 去除 rRNA
        if [ ! -f "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ]; then
            if ls "${RRNA_INDEX}"*.bt2* &> /dev/null; then
                bowtie_tmp_prefix="${BOWTIE_TMP_DIR}/${sample_name}_clean"
                rm -f "${bowtie_tmp_prefix}.1" "${bowtie_tmp_prefix}.2" "${bowtie_tmp_prefix}.1.gz" "${bowtie_tmp_prefix}.2.gz"

                bowtie2 -p ${HIGH_THREADS} --very-fast-local --no-unal -x "${RRNA_INDEX}" \
                        -1 "${TRIM_DIR}/${sample_name}_1.clean.fq.gz" \
                        -2 "${TRIM_DIR}/${sample_name}_2.clean.fq.gz" \
                        --un-conc-gz "${bowtie_tmp_prefix}" \
                        > /dev/null 2> "${CLEAN_DIR}/${sample_name}_bowtie2.log"

                bowtie2_status=$?
                if [ ${bowtie2_status} -ne 0 ]; then
                    echo "❌ [报错] ${sample_name} bowtie2 去 rRNA 失败，详见 ${CLEAN_DIR}/${sample_name}_bowtie2.log"
                    continue
                fi

                # 处理 bowtie2 输出的带后缀文件名
                if [ -f "${bowtie_tmp_prefix}.1" ]; then
                    mv "${bowtie_tmp_prefix}.1" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                elif [ -f "${bowtie_tmp_prefix}.1.gz" ]; then
                    mv "${bowtie_tmp_prefix}.1.gz" "${CLEAN_DIR}/${sample_name}_1.final.fq.gz"
                fi

                if [ -f "${bowtie_tmp_prefix}.2" ]; then
                    mv "${bowtie_tmp_prefix}.2" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
                elif [ -f "${bowtie_tmp_prefix}.2.gz" ]; then
                    mv "${bowtie_tmp_prefix}.2.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz"
                fi
            else
                echo "❌ [报错] rRNA/mtDNA bowtie2 索引未找到: ${RRNA_INDEX}*.bt2*"
                continue
            fi
        fi

        if [ ! -s "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" ] || [ ! -s "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" ]; then
            echo "❌ [报错] ${sample_name} 去 rRNA 后 FASTQ 缺失，跳过 STAR。"
            continue
        fi

        # [3/3] STAR 比对
        if ! valid_bam "${bam_out}"; then
            STAR --runThreadN ${HIGH_THREADS} --genomeDir "${STAR_INDEX}" \
                 --readFilesIn "${CLEAN_DIR}/${sample_name}_1.final.fq.gz" "${CLEAN_DIR}/${sample_name}_2.final.fq.gz" \
                 --readFilesCommand zcat \
                 --outFileNamePrefix "${ALIGN_DIR}/${sample_name}." \
                 --outSAMtype BAM SortedByCoordinate \
                 --winAnchorMultimapNmax 500 --outFilterMultimapNmax 500 \
                 --outMultimapperOrder Random --quantMode GeneCounts --outSAMattributes All \
                 --genomeSAsparseD 3 \
                 --limitBAMsortRAM ${STAR_RAM} > /dev/null

            if valid_bam "${bam_out}"; then
                samtools index -@ ${LOW_THREADS} "${bam_out}"
            else
                echo "❌ [报错] ${sample_name} STAR 未生成有效 BAM。"
                continue
            fi
        fi
    done
else
    echo "❌ [报错] 未在 ${RAW_DIR} 找到 FASTQ 文件。"
fi


# ==========================================
# Step 4: TEcount 定量（家族水平）
# ==========================================
echo "=== Step 4: TEcount 定量 ==="

shopt -s nullglob
bam_files=("${ALIGN_DIR}"/*.Aligned.sortedByCoord.out.bam)
shopt -u nullglob

if [ ${#bam_files[@]} -gt 0 ]; then
    for bam_file in "${bam_files[@]}"; do
        if ! valid_bam "${bam_file}"; then
            echo "⚠️ [跳过] $(basename "$bam_file") 不是有效 BAM。"
            continue
        fi

        sample_name=$(basename "$bam_file" .Aligned.sortedByCoord.out.bam)
        if [ -f "${COUNTS_DIR}/${sample_name}.cntTable" ]; then
            echo "✅ [跳过] ${sample_name} TEcount 已完成。"
            continue
        fi

        (
            TEcount --sortByPos --format BAM --mode multi \
                    --GTF "${GTF_GENE}" --TE "${GTF_TE}" \
                    --project "${COUNTS_DIR}/${sample_name}" \
                    --stranded reverse -b "${bam_file}" \
            && echo "🎉 [TEcount 完成] ${sample_name}"
        ) &
    done
    wait
else
    echo "⚠️ [警告] 未找到 BAM 文件，跳过 TEcount。"
fi

echo "✅ === Step 1-4 完成 (fastp/bowtie2/STAR/TEcount) ==="
