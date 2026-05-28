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
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/')

pbmc1 = qread('Malignant_RNA_assay.qs')

# （注：Seurat 的 NormalizeData 默认是取自然对数 log1p。为了生成 layer = "data" 以防报错，这里保留该步骤）
#计算logCPM
pbmc1 <- NormalizeData(pbmc1,normalization.method = "LogNormalize", scale.factor = 1000000)

# 1. Seurat v5 提取表达矩阵
# 提取 Seurat 里的 data 矩阵
seurat_data = LayerData(pbmc1, assay = "RNA", layer = "data")

# 修正：将 ln(x+1) 转换为 log2(x+1)，只需除以 ln(2)
# 这样 data_log2 才是 CRDscore 要求的 log2 尺度表达矩阵
data_log2 <- as.data.frame(seurat_data / log(2))

# ==========================================
# 定义 Signature 完整路径，并动态提取名称
# ==========================================
# 在这里输入包含文件名和后缀的完整路径
signature_file_path <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/files/CRC_Proliferation_Invasion_Metastasis_Genes.csv"

# 自动提取文件名：去除路径和 .csv 后缀
signature_name <- sub("\\.csv$", "", basename(signature_file_path))

# 【修改点】：设置文件前缀包含 addmodulescore
signature_file_prefix <- paste0(signature_name, "_addmodulescore")

# 2. 提取基因集 (按照完整路径读取)
CRC_data = read.csv(signature_file_path, header = T, check.names = F)
target_genes = as.character(CRC_data[,1])
target_genes = intersect(target_genes, rownames(pbmc1))

# 3. 计算评分 (使用 AddModuleScore)
# AddModuleScore 接收一个基因列表(list)，计算结果会存入 meta.data
pbmc1 <- AddModuleScore(object = pbmc1, features = list(target_genes), name = "ModuleScore")

# 提取评分：Seurat 会在提供的 name 后面自动加上数字 1 (即 ModuleScore1)
score <- pbmc1@meta.data[, "ModuleScore1", drop = FALSE]
colnames(score) <- "score" # 重命名为 score 以兼容后续合并和绘图逻辑
score$id <- rownames(score) # 为后续 merge 提供 ID 列

meta = pbmc1@meta.data
meta$id = rownames(meta) 

cli = read.csv("GSE132465_Cli.csv", header=T, check.names=F)

# 按照患者 ID 合并细胞元数据与临床信息
rt = merge(meta, cli, by.x= "orig.ident",by.y="Tumor")
# 按照细胞 ID 合并评分结果与上一步的合并表
rt1 = merge(score, rt, by.x = "id", by.y = "id")

rt1$metastasis <- ifelse(
  grepl("M0$", rt1$'TNM stage'),
  "Primary",
  "Metastasis"
)

# 【修改点】：动态保存 CSV 文件名
write.csv(rt1, paste0(signature_file_prefix, '_CRC_results.csv'), row.names= FALSE, quote = FALSE)

# ==============================
# 图 1：针对 Stage (分期)
# ==============================
data_stage = rt1[, c('score', 'Stage')]
colnames(data_stage) = c("score", "Type")

# 按照你的要求进行分期合并
data_stage$Type[data_stage$Type %in% c('IIIA', 'IIIB', 'IIIC')] = 'Ⅲ'

# 设置因子水平（去重并排序）
group_stage = levels(factor(data_stage$Type))
data_stage$Type = factor(data_stage$Type, levels = group_stage)

# 生成比较组列表
my_comparisons_stage = list()
if(length(group_stage) >= 2){
  comp_stage = combn(group_stage, 2)
  for(i in 1:ncol(comp_stage)){ my_comparisons_stage[[i]] <- comp_stage[, i] }
}

# 绘图 1 (标题与 Y 轴标签已修改)
p1 = ggplot(data_stage, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle(paste0(signature_name, "+GSE132465+addmodulescore (Stage)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("AddModuleScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p1)
# 动态保存图片 1
ggsave(paste0(signature_file_prefix, "_Stage_plot.png"), plot = p1, width = 6, height = 5, dpi = 300)


# ==============================
# 图 2：针对 Metastasis (转移情况)
# ==============================
data_meta = rt1[, c('score', 'metastasis')]
colnames(data_meta) = c("score", "Type")

# 设置因子水平
group_meta = levels(factor(data_meta$Type))
data_meta$Type = factor(data_meta$Type, levels = group_meta)

# 生成比较组列表
my_comparisons_meta = list()
if(length(group_meta) >= 2){
  comp_meta = combn(group_meta, 2)
  for(i in 1:ncol(comp_meta)){ my_comparisons_meta[[i]] <- comp_meta[, i] }
}

# 绘图 2 (标题与 Y 轴标签已修改)
p2 = ggplot(data_meta, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_meta, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle(paste0(signature_name, "+GSE132465+addmodulescore (Metastasis)")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("AddModuleScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p2)
# 动态保存图片 2
ggsave(paste0(signature_file_prefix, "_Metastasis_plot.png"), plot = p2, width = 6, height = 5, dpi = 300)