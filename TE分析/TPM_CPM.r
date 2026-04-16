# ==============================================================================
# 综合处理脚本：双输出模式 (修复版 - 解决 Y_RNA 重复问题)
# ==============================================================================
setwd('/mnt/windowsdata/qiuzerui/Phf20-26.1.23/')
library(data.table)
library(dplyr)
library(rtracklayer)
library(GenomicRanges)

# ==============================
# 1. 参数与路径设置
# ==============================
counts_file <- "Phf20_1.23_counts.csv"

# 注释文件路径
gene_gtf_path <- "/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38/gencode.vM38.annotation_PRI.gtf"
te_gtf_path   <- "/mnt/windowsdata/qiuzerui/RNAannotations/annotationMv38/m39_TE.gtf"

# 输出文件名
output_cpm_file <- "Phf20_1.23_CPM.csv"  
output_tpm_file <- "Phf20_1.23_TPM.csv"  
#samplename <- c('shNT_rep1','shNT_rep2','shNT_rep3','shPHF20_rep1','shPHF20_rep2','shPHF20_rep3')

# ==============================
# 2. 读取数据与准备注释 (保留唯一ID)
# ==============================
message(paste0("[", Sys.time(), "] 正在读取 Counts 文件..."))
counts_df <- fread(counts_file)
#colnames(counts_df)[2:ncol(counts_df)] <- samplename

# --- 2.1 加载 GTF 提取 Gene Symbol ---
message(paste0("[", Sys.time(), "] 正在加载 GTF 以匹配 Gene Symbol..."))
gene_gtf <- import(gene_gtf_path)
gene_map <- unique(as.data.frame(mcols(gene_gtf)[, c("gene_id", "gene_name")]))

# --- 2.2 关联 Symbol (关键修改：不覆盖 RepeatID) ---
# 使用 left_join 将 symbol 信息并入
counts_df <- left_join(counts_df, gene_map, by = c("RepeatID" = "gene_id"))

# 创建一个新的 Symbol 列：
# 如果匹配到 gene_name，则用 gene_name；如果是 TE (NA)，则暂时填入 RepeatID
counts_df$Symbol <- ifelse(
  !is.na(counts_df$gene_name), 
  counts_df$gene_name, 
  counts_df$RepeatID
)

# 删除多余的 gene_name 列，保留 RepeatID 和 Symbol
counts_df$gene_name <- NULL

# 调整列顺序：把 Symbol 放到 RepeatID 后面
col_order <- c("RepeatID", "Symbol", setdiff(names(counts_df), c("RepeatID", "Symbol")))
counts_df <- counts_df %>% select(all_of(col_order))

message("   - ID 映射完成 (保留 Unique ID, 新增 Symbol 列)")

# ==============================
# 3. 计算并输出 CPM (独立文件)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算 CPM..."))

df_cpm <- copy(counts_df)

# 识别样本列 (排除 ID 和 Symbol)
sample_cols <- setdiff(names(df_cpm), c("RepeatID", "Symbol"))

# 计算 CPM
for (col in sample_cols) {
  library_size <- sum(df_cpm[[col]], na.rm = TRUE)
  cpm_val <- (df_cpm[[col]] / library_size) * 1e6
  df_cpm[[paste0(col, "_CPM")]] <- round(cpm_val, 2)
}

# --- 输出文件 1: CPM ---
message(paste(">>> 正在导出 CPM 文件:", output_cpm_file))
write.csv(df_cpm, output_cpm_file, row.names = FALSE)
message("   ✅ CPM 文件已保存")

# ==============================
# 4. 计算长度 (Gene + TE)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算基因长度 (基于唯一ID)..."))

# --- 4.1 Gene 长度 (使用 gene_id 聚合，确保唯一性) ---
gene_exons <- gene_gtf[gene_gtf$type == "exon"]
# 这里必须用 gene_id 分组，不能用 symbol (否则 Y_RNA 会合并在一起)
gene_exons_list <- split(gene_exons, mcols(gene_exons)$gene_id)
gene_widths <- sum(width(reduce(gene_exons_list)))
gene_len_df <- data.frame(final_id = names(gene_widths), length = as.numeric(gene_widths))

# --- 4.2 TE 长度 ---
te_gtf <- import(te_gtf_path)
te_mcols <- mcols(te_gtf)
# 构造与 Count 表一致的 TE ID (gene_id:family:class)
te_ids <- paste(te_mcols$gene_id, te_mcols$family_id, te_mcols$class_id, sep = ":")
mcols(te_gtf)$te_unique_id <- te_ids
te_list <- split(te_gtf, mcols(te_gtf)$te_unique_id)
te_widths <- sum(width(reduce(te_list)))
te_len_final <- data.frame(final_id = names(te_widths), length = as.numeric(te_widths))

# --- 4.3 合并长度表 ---
all_lengths <- rbind(gene_len_df, te_len_final)

# ==============================
# 5. 计算并输出 TPM (独立文件)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算 TPM..."))

# --- 5.1 合并长度 (基于 RepeatID <-> final_id) ---
# 此时 counts_df 里的 RepeatID 是唯一的 ENSG 或 TE_ID，完全匹配 length 表
df_tpm <- left_join(counts_df, all_lengths, by = c("RepeatID" = "final_id"))

# 移除无长度的行
df_tpm <- df_tpm[!is.na(df_tpm$length), ]

# 计算 TPM 函数
calculate_tpm <- function(counts, lengths) {
  # 防止长度为0导致的错误
  lengths <- ifelse(lengths == 0, 1, lengths)
  rpk <- counts / (lengths / 1000)
  # 避免 sum(rpk) 为 0
  total_rpk <- sum(rpk, na.rm = TRUE)
  scaling_factor <- ifelse(total_rpk == 0, 1, total_rpk / 1e6)
  return(rpk / scaling_factor)
}

# 仅对样本列计算
tpm_calc_cols <- setdiff(names(df_tpm), c("RepeatID", "Symbol", "length"))

for (col in tpm_calc_cols) {
  tpm_val <- calculate_tpm(df_tpm[[col]], df_tpm$length)
  df_tpm[[paste0(col, "_TPM")]] <- round(tpm_val, 4)
}

# --- 5.2 整理输出列顺序 ---
# 顺序：RepeatID, Symbol, length, Raw_Counts..., TPM_Values...
cols_meta <- c("RepeatID", "Symbol", "length")
# 如果只想保留 TPM 值，把下面一行改为只选包含 "_TPM" 的列
cols_data <- setdiff(names(df_tpm), cols_meta) 

df_tpm_out <- df_tpm %>% select(all_of(c(cols_meta, cols_data)))

# --- 输出文件 2: TPM ---
message(paste(">>> 正在导出 TPM 文件:", output_tpm_file))
write.csv(df_tpm_out, output_tpm_file, row.names = FALSE)
message("   ✅ TPM 文件已保存")

message("========================================================")
message("🎉 全部完成！已生成文件：")
message(paste("1.", output_cpm_file))
message(paste("2.", output_tpm_file))
message("========================================================")
# ==============================================================================
# TE 差异分析脚本 (最终修复版)
# 功能：
# 1. 筛选 TE
# 2. 拆分 SubFamily, Family, Class
# 3. 计算 Log2FC 和 PValue
# 4. 仅保留分类、统计值和 CPM 数据 (去除原始 RepeatID, Symbol, Counts)
# ==============================================================================

# 加载必要的包
library(data.table)
library(dplyr)
library(tidyr)
library(stringr) # 确保加载 stringr

# --- 设置输入输出文件 ---
input_cpm <- "Phf20_1.23_CPM.csv"
output_te_stats <- "Phf20_1.23_TE_Stats_Clean.csv"

message(paste0("[", Sys.time(), "] 正在读取 CPM 文件..."))
df <- fread(input_cpm)

# 1. 筛选 TE (至少包含2个冒号)
#    确保使用 stringr 包的函数
df_te <- df[stringr::str_count(df$RepeatID, ":") >= 2, ]
message(paste0("   - 筛选出 TE 数量: ", nrow(df_te)))

# 2. 拆分 ID 为三列 (SubFamily, Family, Class)
#    注意：separate 后转为 data.frame 避免 data.table 的索引兼容性问题
df_te <- df_te %>%
  separate(RepeatID, c("SubFamily", "Family", "Class"), sep = ":", remove = FALSE, extra = "merge") %>%
  as.data.frame()

# 3. 定义 CPM 列名 (用于计算和输出)
ctrl_cpm_cols <- c("L1MKL2609676-Scr_1_Mixt_CPM", "L1MKL2609677-Scr_2_Mixt_CPM")
treat_cpm_cols <- c("L1MKL2609678-Phf20_1_Mixt_CPM", "L1MKL2609679-Phf20_2_Mixt_CPM")
all_cpm_cols <- c(ctrl_cpm_cols, treat_cpm_cols)

# 检查列名是否存在
if (!all(all_cpm_cols %in% names(df_te))) {
  stop("错误：CPM 列名匹配失败，请检查文件中的列名！")
}

# ==============================
# 4. 计算 Log2FC (向量化计算)
# ==============================
message(paste0("[", Sys.time(), "] 正在计算统计指标..."))

# 提取 CPM 数据矩阵
mat_ctrl <- df_te[, ctrl_cpm_cols]
mat_treat <- df_te[, treat_cpm_cols]

# 计算均值 (添加 0.01 防止 log(0))
pseudo <- 0.01
mean_ctrl <- rowMeans(mat_ctrl, na.rm = TRUE)
mean_treat <- rowMeans(mat_treat, na.rm = TRUE)

# 计算 Log2FC
df_te$Log2FC <- round(log2(mean_treat + pseudo) - log2(mean_ctrl + pseudo), 4)

# ==============================
# 5. 计算 P-value (定义函数并计算)
# ==============================

# 定义 T-test 函数 (确保这是一个完整的代码块)
calc_pval <- function(x, c_cols, t_cols) {
  v_c <- as.numeric(x[c_cols])
  v_t <- as.numeric(x[t_cols])
  
  # 如果方差均为0 (数值完全一致)，直接返回 1
  if (var(v_c) == 0 && var(v_t) == 0) return(1)
  
  tryCatch({
    t.test(v_t, v_c)$p.value
  }, error = function(e) NA)
}

# 执行计算 (apply 会逐行调用上述函数)
df_te$PValue <- apply(df_te, 1, calc_pval, c_cols = ctrl_cpm_cols, t_cols = treat_cpm_cols)
df_te$PValue <- round(df_te$PValue, 5)

# ==============================
# 6. 整理输出 (移除 RepeatID, Symbol, Counts)
# ==============================

# 定义最终要保留的列：分类 + 统计 + CPM
final_cols <- c(
  "SubFamily", "Family", "Class", 
  "Log2FC", "PValue", 
  all_cpm_cols
)

# 筛选列
df_final <- df_te %>% select(all_of(final_cols))

# 导出文件
message(paste(">>> 正在导出精简版统计文件:", output_te_stats))
write.csv(df_final, output_te_stats, row.names = FALSE)

message("   ✅ 全部完成！已生成文件：", output_te_stats)