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

# --- 读取并构建 Seurat 对象 ---


# 替换 UMI 矩阵读取
# data.table 读取后默认是 data.table 格式，我们需要转回 data.frame 并设置行名
counts_matrix <- fread("GSE132465_10X_UMI_matrix.txt", header = TRUE, sep = "\t", check.names = FALSE)
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- counts_matrix[[1]] # 第一列设为行名
counts_matrix <- counts_matrix[,-1]          # 删掉原来的第一列

# 替换注释文件读取
cell_annotation <- fread("GSE132465_10X_cell_annotation.txt", header = TRUE, sep = "\t")
cell_annotation <- as.data.frame(cell_annotation)
rownames(cell_annotation) <- cell_annotation[[1]]
cell_annotation <- cell_annotation[,-1]

scRNA <- CreateSeuratObject(counts = counts_matrix, meta.data = cell_annotation)
pbmc1 <- subset(scRNA, subset = Class == 'Tumor')

# （注：Seurat 的 NormalizeData 默认是取自然对数 log1p。为了生成 layer = "data" 以防报错，这里保留该步骤）
pbmc1 <- NormalizeData(pbmc1)

# 1. Seurat v5 提取表达矩阵
# 提取 Seurat 里的 data 矩阵
seurat_data = LayerData(pbmc1, assay = "RNA", layer = "data")

# 修正：将 ln(x+1) 转换为 log2(x+1)，只需除以 ln(2)
# 这样 data_log2 才是 CRDscore 要求的 log2 尺度表达矩阵
data_log2 <- as.data.frame(seurat_data / log(2))

# 2. 提取基因集
CRC_data = read.csv("/mnt/disk1/qiuzerui/downloads/CRC/signature/112 primary cilium genes.csv", header = T, check.names = F)

target_genes = as.character(CRC_data[,1])
target_genes = intersect(target_genes, rownames(pbmc1))

# 3. 计算评分
# 【后续用 data_log2 替换】：计算得分时的 expr 参数改为 data_log2
score <- cal_CRDscore(expr = data_log2, n.bins = 50, circadians = 
                        target_genes, study.type = "scRNAseq")
gc()
score = as.data.frame(score)
score = cbind(id=rownames(score), score)

score$id = rownames(score)
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
write.csv(rt1,'112 primary cilium genes_CRC_CRDscore.csv',row.names= FALSE,quote = FALSE)
# ==============================
# 图 1：针对 Stage (分期)
# ==============================
data_stage = rt1[, c('score', 'Stage')]
colnames(data_stage) = c("score", "Type")

# 按照你的要求进行分期合并
data_stage$Type[data_stage$Type %in% c('IIIA', 'IIIB', 'IIIC')] = 'Ⅲ'
data_stage$Type[data_stage$Type == 'IIA'] = 'I'
# 如果有其他分期（如 II），建议也确认一下是否需要保留

# 设置因子水平（去重并排序）
group_stage = levels(factor(data_stage$Type))
data_stage$Type = factor(data_stage$Type, levels = group_stage)

# 生成比较组列表
my_comparisons_stage = list()
if(length(group_stage) >= 2){
  comp_stage = combn(group_stage, 2)
  for(i in 1:ncol(comp_stage)){ my_comparisons_stage[[i]] <- comp_stage[, i] }
}

# 绘图 1
p1 = ggplot(data_stage, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_stage, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle("112 primary cilium genes+GSE123465+CRDscore (Stage)") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p1)
ggsave("112_primary_cilium_genes_Stage_plot.png", plot = p1, width = 6, height = 5, dpi = 300)


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

# 绘图 2
p2 = ggplot(data_meta, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom = "errorbar", width = 0.6) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  stat_compare_means(comparisons = my_comparisons_meta, method = "wilcox.test") +
  theme(legend.position = "none") +
  ggtitle("112 primary cilium genes+GSE123465+CRDscore (Metastasis)") +
  coord_cartesian(ylim = c(-0.2, 0.2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab("CRDScore") +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.text.x = element_text(size = 15, face = "bold", color = "black"),
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p2)
ggsave("112_primary_cilium_genes_Metastasis_plot.png", plot = p2, width = 6, height = 5, dpi = 300)