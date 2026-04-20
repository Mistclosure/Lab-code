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

# --- 承接预处理脚本：设置工作目录并读取 .qs 对象 ---
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE231559')

# 1. 读取上一步预处理完成的高质量 Seurat 对象
pbmc1 <- qread("Malignant_RNA_assay.qs")

# 2. 筛选肿瘤细胞（注意：预处理脚本中标记肿瘤的列名为 Type，取值为 Tumor）

# （Seurat NormalizeData 默认 log1p 处理，生成 layer = "data"）
pbmc1 <- NormalizeData(pbmc1)

# ------------------------------------------------------------------------------
# 3. 提取表达矩阵并进行 log2 转换
# ------------------------------------------------------------------------------
# Seurat v5 提取 data 层矩阵
seurat_data = LayerData(pbmc1, assay = "RNA", layer = "data")

# 修正：将 ln(x+1) 转换为 log2(x+1)，只需除以 ln(2)
# 这样 data_log2 才是 CRDscore 要求的 log2 尺度表达矩阵
data_log2 <- as.data.frame(seurat_data / log(2))
#仅第一次运行时使用，用于储存结果文件
pbmc1[['RNA']]$data_log2 = data_log2
qsave(pbmc1,'Malignant_GSE231559.qs')
# ------------------------------------------------------------------------------
# 4. 计算评分
# ------------------------------------------------------------------------------
# --- 动态提取基因集文件名 ---
signature_file <- "/mnt/disk1/qiuzerui/downloads/CRC/signature/112 primary cilium genes.csv"
signature_name <- tools::file_path_sans_ext(basename(signature_file))

CRC_data = read.csv(signature_file, header = T, check.names = F)

target_genes = as.character(CRC_data[,1])
target_genes = intersect(target_genes, rownames(pbmc1))

# 计算得分，使用准备好的 data_log2
score <- cal_CRDscore(expr = data_log2, n.bins = 50, circadians = 
                        target_genes, study.type = "scRNAseq")
gc()

score = as.data.frame(score)
score$id = rownames(score)

# ------------------------------------------------------------------------------
# 5. 合并元数据与临床信息
# ------------------------------------------------------------------------------
meta = pbmc1@meta.data
meta$id = rownames(meta) 

# 加载临床信息
cli = read.csv("GSE231559_Cli.csv", header=T, check.names=F) 

# 按照 Patient ID 或 orig.ident 合并
rt = merge(meta, cli,  by="Patient") 

# 按照细胞 ID 合并评分结果
rt1 = merge(score, rt, by.x = "id", by.y = "id")

# 动态保存 CSV 文件名
write.csv(rt1, paste0(signature_name, '_CRC_CRDscore.csv'), row.names= FALSE, quote = FALSE)

# ==============================================================================
# 图 1：针对 Stage (分期)
# ==============================================================================
data_stage = rt1[, c('score', 'Stage')]
colnames(data_stage) = c("score", "Type")

# 分期合并逻辑

group_stage = levels(factor(data_stage$Type))
data_stage$Type = factor(data_stage$Type, levels = group_stage)

my_comparisons_stage = list()
if(length(group_stage) >= 2){
  comp_stage = combn(group_stage, 2)
  for(i in 1:ncol(comp_stage)){ my_comparisons_stage[[i]] <- comp_stage[, i] }
}

p1 = ggplot(data_stage, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test") +
  theme(legend.position = "none") +
  # 动态标题
  ggtitle(paste0(signature_name, "+GSE231559+CRDscore (Stage)")) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p1)
# 动态保存图片名
ggsave(paste0(signature_name, "_Stage_plot.png"), plot = p1, width = 6, height = 5, dpi = 300)

# ==============================================================================
# 图 2：针对 Metastasis (转移情况)
# ==============================================================================
data_meta = rt1[, c('score', 'metastasis')]
colnames(data_meta) = c("score", "Type")

group_meta = levels(factor(data_meta$Type))
data_meta$Type = factor(data_meta$Type, levels = group_meta)

my_comparisons_meta = list()
if(length(group_meta) >= 2){
  comp_meta = combn(group_meta, 2)
  for(i in 1:ncol(comp_meta)){ my_comparisons_meta[[i]] <- comp_meta[, i] }
}

p2 = ggplot(data_meta, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_meta, method = "wilcox.test") +
  theme(legend.position = "none") +
  # 动态标题
  ggtitle(paste0(signature_name, "+GSE231559+CRDscore (Metastasis)")) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p2)
# 动态保存图片名
ggsave(paste0(signature_name, "_Metastasis_plot.png"), plot = p2, width = 6, height = 5, dpi = 300)
