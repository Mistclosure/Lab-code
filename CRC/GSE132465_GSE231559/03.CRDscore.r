# ==============================================================================
# 单细胞 CRDscore 评分与临床相关性分析 (适配整合后的 Merged_malig.qs)
# ==============================================================================

# 1. 载入所需的 R 包 -----------------------------------------------------------
library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(tidyverse)
library(plyr)
library(Seurat)
library(CRDscore)
library(qs)
library(data.table)

# 设置工作目录至整合数据集所在路径
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465_GSE231559/')

# 2. 读取数据与预处理 ----------------------------------------------------------
cat("正在读取合并后的恶性细胞子集...\n")
malig_seurat <- qread('Merged_malig.qs')

# 直接从 Seurat v5 中提取 data 图层
seurat_data <- LayerData(malig_seurat, assay = "RNA", layer = "data")

# 转换为 log2 尺度表达矩阵 (通过除以 ln(2) 实现 CRDscore 要求的 log2 尺度)
expr_log2 <- as.data.frame(seurat_data / log(2))

# ==============================================================================
# 3. 提取基因集并计算 CRDscore -------------------------------------------------
# ==============================================================================
signature_name <- "ciliopathy_genes"
signature_file_prefix <- signature_name # 若不需要空格可改用 gsub(" ", "_", signature_name)

cat("正在读取 Signature 基因集...\n")
CRC_data <- read.csv(paste0("/mnt/disk1/qiuzerui/downloads/CRC/signature/", signature_name, ".csv"), header = TRUE, check.names = FALSE)
target_genes <- as.character(CRC_data[, 1])
target_genes <- intersect(target_genes, rownames(malig_seurat))

cat("正在计算 CRDscore...\n")
crd_score <- cal_CRDscore(expr = expr_log2, n.bins = 50, circadians = target_genes, study.type = "scRNAseq")
gc()

# 整理评分数据框
crd_score_df <- as.data.frame(crd_score)
crd_score_df <- cbind(id = rownames(crd_score_df), crd_score_df)
colnames(crd_score_df)[2] <- "score"
crd_score_df$id <- rownames(crd_score_df)

# ==============================================================================
# 4. 临床信息合并 --------------------------------------------------------------
# ==============================================================================
cat("正在合并临床信息...\n")
meta_data <- malig_seurat@meta.data
meta_data$id <- rownames(meta_data) 

# 直接将评分结果与包含完整临床信息的 meta_data 合并
merged_score_meta <- merge(crd_score_df, meta_data, by = "id")

# 动态保存 CSV 文件
write.csv(merged_score_meta, paste0(signature_file_prefix, '_Merged_CRC_CRDscore.csv'), row.names = FALSE, quote = FALSE)

# ==============================================================================
# 5. 可视化分析 ----------------------------------------------------------------
# ==============================================================================
cat("正在生成并保存箱线图...\n")

# ------------------------------
# 图 1：针对 Stage (分期)
# ------------------------------
# 直接使用已有的 Stage 列
data_stage <- merged_score_meta[, c('score', 'Stage')]
colnames(data_stage) <- c("score", "Type")

group_stage <- levels(factor(data_stage$Type))
data_stage$Type <- factor(data_stage$Type, levels = group_stage)

my_comparisons_stage <- list()
if(length(group_stage) >= 2){
  comp_stage <- combn(group_stage, 2)
  for(i in 1:ncol(comp_stage)){ my_comparisons_stage[[i]] <- comp_stage[, i] }
}

p1 <- ggplot(data_stage, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle(paste0(signature_name, " Merged_CRC CRDscore (Stage)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p1)
ggsave(paste0(signature_file_prefix, "_Stage_plot.png"), plot = p1, width = 6, height = 5, dpi = 300)

# ------------------------------
# 图 2：针对 Metastasis (转移情况)
# ------------------------------
# 直接使用已有的 metastasis 列，去除了多余的 ifelse 判断
data_meta <- merged_score_meta[, c('score', 'metastasis')]
colnames(data_meta) <- c("score", "Type")

group_meta <- levels(factor(data_meta$Type))
data_meta$Type <- factor(data_meta$Type, levels = group_meta)

my_comparisons_meta <- list()
if(length(group_meta) >= 2){
  comp_meta <- combn(group_meta, 2)
  for(i in 1:ncol(comp_meta)){ my_comparisons_meta[[i]] <- comp_meta[, i] }
}

p2 <- ggplot(data_meta, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_meta, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle(paste0(signature_name, " Merged_CRC CRDscore (Metastasis)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p2)
ggsave(paste0(signature_file_prefix, "_Metastasis_plot.png"), plot = p2, width = 6, height = 5, dpi = 300)

cat("全部分析并绘图完毕！\n")