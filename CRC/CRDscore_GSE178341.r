library(ggplot2)
library(ggpubr)
library(limma)
library(reshape2)
library(tidyverse)
library(plyr)
library(Seurat)
library(CRDscore)
library(qs)

# --- 基础配置 ---
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE178341/')
dir.create("pictures", showWarnings = FALSE) 

# 定义输入文件名变量DNA-damage-response genes.csv   ciliopathy_genes.csv
input_file <- "/mnt/disk1/qiuzerui/downloads/CRC/signature/112 primary cilium genes.csv"
file_base_name <- tools::file_path_sans_ext(basename(input_file))

pbmc1 = qread('Malignant.qs')
# 过滤全0基因
pbmc1 <- pbmc1[Matrix::rowSums(GetAssayData(pbmc1, layer = "counts")) > 0, ]
# 1. Seurat v5 提取表达矩阵
exp = as.data.frame(LayerData(pbmc1, assay = "RNA", layer = "data"))

# 2. 提取基因集
CRC_data = read.csv(input_file, header = T, check.names = F)
target_genes = as.character(CRC_data[,1])
target_genes = intersect(target_genes, rownames(pbmc1))

# --- 3. 计算评分 (方法 A: CRDscore) ---
score <- cal_CRDscore(expr = exp, n.bins = 50, circadians = 
                        target_genes, study.type = "scRNAseq")
gc()
score = as.data.frame(score)
score$id = rownames(score)

# --- 4. 计算评分 (方法 B: AddModuleScore) ---
# Seurat 会将结果存入 meta.data，列名为 name + 1
pbmc1 <- AddModuleScore(pbmc1, features = list(target_genes), name = "AddModuleScore")

# 提取 metadata (包含 AddModuleScore1)
meta = pbmc1@meta.data
meta$id = rownames(meta) 
cli = read.csv("Cli.csv", header=T, check.names=F)

# 1. 按照患者 ID 合并临床信息
rt = merge(meta, cli, by.x= "orig.ident",by.y="PatientBarcode")

# 2. 按照细胞 ID 合并所有结果
rt1 = merge(score, rt, by.x = "id", by.y = "id")

# --- 绘图公共参数设置 ---
type_group = levels(factor(rt1$"Tumor Stage"))
my_comparisons = list()
if(length(type_group) >= 2){
  comp = combn(type_group, 2)
  for(i in 1:ncol(comp)){ my_comparisons[[i]] <- comp[,i] }
}
rt1$'Tumor Stage'[rt1$'Tumor Stage'==2] =1
# ==================== 绘图 1: CRDscore ====================
data_p = rt1[,c('score','Tumor Stage')]
colnames(data_p) = c("score","Type")
data_p$Type = factor(data_p$Type, levels=type_group)

p = ggplot(data_p, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom="errorbar", width=0.6) +
  geom_boxplot(alpha=0.7, outlier.shape = NA, size=0.7, width=0.7, fatten=0.7) +
  theme_bw() + theme(panel.grid=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ggtitle(paste0(file_base_name, "+CRC+CRDscore")) +
  ylim(-0.2, 0.2) + ylab("CRDScore") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p)
ggsave(paste0("pictures/", file_base_name, "_CRC_CRDscore.png"), plot = p, width = 7, height = 6, dpi = 300)

# ==================== 绘图 2: AddModuleScore ====================
# 注意：AddModuleScore 的列名通常是 "AddModuleScore1"
data_p2 = rt1[,c('AddModuleScore1','Tumor Stage')]
colnames(data_p2) = c("score","Type")
data_p2$Type = factor(data_p2$Type, levels=type_group)

p2 = ggplot(data_p2, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom="errorbar", width=0.6) +
  geom_boxplot(alpha=0.7, outlier.shape = NA, size=0.7, width=0.7, fatten=0.7) +
  theme_bw() + theme(panel.grid=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ggtitle(paste0(file_base_name, "+CRC+AddModuleScore")) +
  # AddModuleScore 的数值范围通常比 CRDscore 大，这里取消 ylim 或根据实际调整
  ylab("AddModuleScore") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p2)
ggsave(paste0("pictures/", file_base_name, "_CRC_AddModuleScore.png"), plot = p2, width = 7, height = 6, dpi = 300)
rownames(rt1) <- rt1$id
pbmc1$type <- rt1[colnames(pbmc1), "Liver metastasis (LM)"]
VlnPlot(pbmc1, features = c("TP53", "BRCA1"), group.by = "type")
# 选出你 112 个基因中在 pbmc1 里真正存在的那些
genes_to_plot <- intersect(target_genes, rownames(pbmc1))

# 画气泡图：颜色代表表达量高低，点的大小代表表达细胞比例
# 1. 快速找两组间的差异基因（仅针对这 112 个基因）
de_results <- FindMarkers(pbmc1, 
                          ident.1 = "primary", 
                          ident.2 = "metastasis", 
                          group.by = "type",
                          features = target_genes,
                          logfc.threshold = 0)

# 2. 按照 Log2FC 的绝对值排序，选出前 30 个
top_30_genes <- rownames(de_results[order(abs(de_results$avg_log2FC), decreasing = TRUE), ])[1:30]

# 3. 只画这 30 个基因
p3=DotPlot(pbmc1, features = top_30_genes, group.by = "type") + coord_flip()
ggsave(paste0("pictures/", file_base_name, "_CRC_dimplot.png"), plot = p3, width = 7, height = 20, dpi = 300,bg = "white")
