#!/bin/bash

# =====================================================================
# 路径与配置区域
# =====================================================================
# 1. 原始数据源目录（支持存放 .sra 实体或已经解压好的华大 4 文件）
SRC_DIR="/mnt/disk1/qiuzerui/expriments/cortex_Y90C/rawdata"

# 2. 工作与输出总目录
OUT_DIR="/mnt/disk1/qiuzerui/expriments/cortex_Y90C"

# 3. 规范化：单独创建测序数据隔离存放目录（存放标准化命名的软链接或解压文件）
FASTQ_DIR="$OUT_DIR/fastqs"

# 4. 参考基因组路径
# ⚠️【核心警告】：此路径必须是通过华大 v3.0 的 `dnbc4tools rna mkref` 建立的专用索引目录！
REF_PATH="/mnt/windowsdata/qiuzerui/scannotations/mouse/mgi_mm10_v3_index/mm10"

# 5. Apptainer 容器沙盒镜像路径
SINGULARITY_IMAGE="/mnt/disk1/qiuzerui/apptainer/dnbc4tools.sif"


# =====================================================================
# 环境变量与资源配置 (高并发生产模式)
# =====================================================================
# 资源分配策略
THREADS=150          # dnbc4tools 运行的最大核心数
DUMP_THREADS=32      # fasterq-dump 线程限制在 32，防止硬盘 I/O 瞬间阻塞锁死
ZIP_THREADS=32       # pigz 后台并发压缩线程


# =====================================================================
# 自动化样本名检测与隔离目录构建 (高鲁棒性清洗算法)
# =====================================================================
mkdir -p "$FASTQ_DIR"
cd "$OUT_DIR"

echo "========================================================"
echo "当前工作主目录: $(pwd)"
echo "标准化 FASTQ 存储目录: $FASTQ_DIR"
echo "正在扫描原始数据目录，自动提取并构建唯一样本列表..."

SAMPLES_RAW=()

if [ ! -d "$SRC_DIR" ]; then
    echo "【错误】原始数据源目录 $SRC_DIR 不存在，请检查路径！"
    exit 1
fi

# 升级版核心过滤逻辑：支持减号(-)与下划线(_)混合型复杂的华大单细胞命名清洗
for file in "$SRC_DIR"/*; do
    if [ -f "$file" ]; then
        filename=$(basename "$file")
        
        if [[ "$filename" == *.sra ]]; then
            # .sra 文件直接提取前缀
            SAMPLES_RAW+=("${filename%.sra}")
            
        elif [[ "$filename" == *.fastq.gz || "$filename" == *.fq.gz ]]; then
            # 剥离压缩后缀
            base="${filename%.fastq.gz}"
            base="${base%.fq.gz}"
            
            # 【精准修正】：将 [R1-4] 修正为 R[1-4]，防止误切样本名自带的 -1 和 -2 
            # 示例：kat8-P60-WT-1-cDNA_R1 -> kat8-P60-WT-1
            sample_name=$(echo "$base" | sed -E 's/[-_](cDNA|Oligo|cdna|oligo).*//; s/[-_]R[1-4].*//; s/_[1-4].*//')
            
            if [ -n "$sample_name" ]; then
                SAMPLES_RAW+=("$sample_name")
            fi
        fi
    fi
done

# 对提取出来的样本列表进行去重和排序
if [ ${#SAMPLES_RAW[@]} -eq 0 ]; then
    echo "【错误】未在 $SRC_DIR 中扫描到任何有效的 .sra 或 .fastq.gz 文件！"
    exit 1
fi
SAMPLES=($(printf "%s\n" "${SAMPLES_RAW[@]}" | sort -u))

echo "成功识别并加载以下 ${#SAMPLES[@]} 个核心单细胞样本:"
for s in "${SAMPLES[@]}"; do echo "  ➔  $s"; done
echo "========================================================"


# =====================================================================
# 主循环程序 (兼容 SRA 4文件切分 与 隔离保护双轨重组)
# =====================================================================
echo "开始串行处理多样本任务..."

for SAMPLE in "${SAMPLES[@]}"; do
    echo "--------------------------------------------------------"
    echo ">>>> 正在处理样本: $SAMPLE"
    
    # 定义当前样本的独立定量结果输出路径
    SAMPLE_OUT_DIR="$OUT_DIR/Output_${SAMPLE}"
    
    # --- 步骤 1: 测序数据 4 文件标准化整合 ---
    if [ ! -f "$FASTQ_DIR/${SAMPLE}_cDNA_R1.fastq.gz" ] || \
       [ ! -f "$FASTQ_DIR/${SAMPLE}_cDNA_R2.fastq.gz" ] || \
       [ ! -f "$FASTQ_DIR/${SAMPLE}_Oligo_R1.fastq.gz" ] || \
       [ ! -f "$FASTQ_DIR/${SAMPLE}_Oligo_R2.fastq.gz" ]; then
        
        # 情况 A: 目录下存在对应的 .sra 实体文件
        if [ -f "$SRC_DIR/${SAMPLE}.sra" ]; then
            echo "[$SAMPLE] 状态: 检测到 SRA 源文件，启动 4 文件高速拆解..."
            
            fasterq-dump --split-files --include-technical \
                         -e $DUMP_THREADS \
                         -O "$FASTQ_DIR" \
                         -t "$FASTQ_DIR/tmp_$SAMPLE" \
                         "$SRC_DIR/${SAMPLE}.sra"
            
            echo "[$SAMPLE] 正在启动高速多线程后台并发压缩 (pigz)..."
            pigz -p $ZIP_THREADS "$FASTQ_DIR/${SAMPLE}_1.fastq" && mv "$FASTQ_DIR/${SAMPLE}_1.fastq.gz" "$FASTQ_DIR/${SAMPLE}_cDNA_R1.fastq.gz" &
            pigz -p $ZIP_THREADS "$FASTQ_DIR/${SAMPLE}_2.fastq" && mv "$FASTQ_DIR/${SAMPLE}_2.fastq.gz" "$FASTQ_DIR/${SAMPLE}_cDNA_R2.fastq.gz" &
            pigz -p $ZIP_THREADS "$FASTQ_DIR/${SAMPLE}_3.fastq" && mv "$FASTQ_DIR/${SAMPLE}_3.fastq.gz" "$FASTQ_DIR/${SAMPLE}_Oligo_R1.fastq.gz" &
            pigz -p $ZIP_THREADS "$FASTQ_DIR/${SAMPLE}_4.fastq" && mv "$FASTQ_DIR/${SAMPLE}_4.fastq.gz" "$FASTQ_DIR/${SAMPLE}_Oligo_R2.fastq.gz" &
            wait
            rm -rf "$FASTQ_DIR/tmp_$SAMPLE"

        # 情况 B: 目录下是已经解压好的原始 4 个 fastq.gz/fq.gz 文件
        else
            echo "[$SAMPLE] 状态: 检测到源目录为 Fastq 压缩包，正在智能检索读段分配..."
            
            # 抓取该核心样本的所有压缩包并排序
            sample_files=($(ls "$SRC_DIR"/${SAMPLE}* 2>/dev/null | grep -E "\.(fastq\.gz|fq\.gz)$" | sort))
            
            if [ ${#sample_files[@]} -ne 4 ]; then
                echo "【严重错误】未能在 $SRC_DIR 中找齐 $SAMPLE 对应的 4 个配对数据！实际仅找到 ${#sample_files[@]} 个。"
                echo "找齐的关联文件列表如下:"
                for f_err in "${sample_files[@]}"; do echo "  -> $(basename $f_err)"; done
                exit 1
            fi

            src_c1="" src_c2="" src_o1="" src_o2=""
            SAMPLE_LOWER=$(echo "$SAMPLE" | tr '[:upper:]' '[:lower:]')
            
            # 前缀剥离算法，完美避免样本名自带的数字（如 WT-1）导致读段识别错乱
            for f in "${sample_files[@]}"; do
                f_base=$(basename "$f")
                f_lower=$(echo "$f_base" | tr '[:upper:]' '[:lower:]')
                
                # 核心安全操作：只截取样本名之后的文件尾缀进行特征识别
                tail_part="${f_lower#*${SAMPLE_LOWER}}"
                
                if [[ "$tail_part" == *"cdna"* ]]; then
                    if [[ "$tail_part" =~ (r1|[-_]1) ]]; then src_c1="$f"; fi
                    if [[ "$tail_part" =~ (r2|[-_]2) ]]; then src_c2="$f"; fi
                elif [[ "$tail_part" == *"oligo"* ]]; then
                    if [[ "$tail_part" =~ (r1|[-_]1) ]]; then src_o1="$f"; fi
                    if [[ "$tail_part" =~ (r2|[-_]2) ]]; then src_o2="$f"; fi
                fi
            done

            # 安全退化机制：如果文件名极其抽象不含任何显式字母，则按字母默认正序进行强制映射
            if [ -z "$src_c1" ] || [ -z "$src_c2" ] || [ -z "$src_o1" ] || [ -z "$src_o2" ]; then
                echo "[$SAMPLE] 提示: 尾缀不规则，已自动启动物理位置顺序重组..."
                src_c1="${sample_files[0]}"
                src_c2="${sample_files[1]}"
                src_o1="${sample_files[2]}"
                src_o2="${sample_files[3]}"
            fi

            # 跨目录建立标准化命名软链接（满足规范且不吃双倍磁盘空间）
            ln -sf "$src_c1" "$FASTQ_DIR/${SAMPLE}_cDNA_R1.fastq.gz"
            ln -sf "$src_c2" "$FASTQ_DIR/${SAMPLE}_cDNA_R2.fastq.gz"
            ln -sf "$src_o1" "$FASTQ_DIR/${SAMPLE}_Oligo_R1.fastq.gz"
            ln -sf "$src_o2" "$FASTQ_DIR/${SAMPLE}_Oligo_R2.fastq.gz"
            echo "[$SAMPLE] 标准化 4 文件隔离软链接构建成功。"
        fi
    else
        echo "[$SAMPLE] 标准化 4 文件 FASTQ 数据已就绪，跳过整合阶段。"
    fi

    # --- 步骤 2: 华大 v3.0 dnbc4tools 定量分析 (Apptainer 容器穿透模式) ---
    if [ ! -d "$SAMPLE_OUT_DIR" ]; then
        echo "[$SAMPLE] 步骤 2: 启动华大 dnbc4tools rna run v3.0 定量流程 (容器版)..."
        
        # 穿透运行：调用容器内部物理锁死的分析引擎，全面映射宿主机路径
        apptainer exec --bind /mnt:/mnt "$SINGULARITY_IMAGE" dnbc4tools rna run \
            --name "$SAMPLE" \
            --cDNAfastq1 "$FASTQ_DIR/${SAMPLE}_cDNA_R1.fastq.gz" \
            --cDNAfastq2 "$FASTQ_DIR/${SAMPLE}_cDNA_R2.fastq.gz" \
            --oligofastq1 "$FASTQ_DIR/${SAMPLE}_Oligo_R1.fastq.gz" \
            --oligofastq2 "$FASTQ_DIR/${SAMPLE}_Oligo_R2.fastq.gz" \
            --genomeDir "$REF_PATH" \
            --threads $THREADS \
            --outdir "$SAMPLE_OUT_DIR"
    else
        echo "[$SAMPLE] Output 定量目录已存在，跳过该样本的定量。"
    fi
    
    echo "<<<< 样本 $SAMPLE 华大自动化流程安全结束！"
done

echo "--------------------------------------------------------"
echo "🎉 所有样本全自动化处理完成！最终规整成果位于各个 Output 目录中。"