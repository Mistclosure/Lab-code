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
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE178318/')
PLOTS_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/cholesterol_Score/GSE178318"
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

file_base_name <- "cholesterolscore"

pbmc1 = qread('Malignant.qs')
pbmc1 <- NormalizeData(pbmc1,normalization.method = "LogNormalize", scale.factor = 1000000)
pbmc1$type <- sub(".*_", "", colnames(pbmc1))
pbmc1 = subset(pbmc1,subset = pbmc1$type %in% c('CRC','LM'))
#pbmc1 = subset(pbmc1,subset = pbmc1$orig.ident == 'COL17')
# 过滤全0基因
pbmc1 <- pbmc1[Matrix::rowSums(GetAssayData(pbmc1, layer = "counts")) > 0, ]
# 1. Seurat v5 提取表达矩阵
# 提取 Seurat 里的 data 矩阵
seurat_data = LayerData(pbmc1, assay = "RNA", layer = "data")

# 修正：将 ln(x+1) 转换为 log2(x+1)，只需除以 ln(2)
# 这样 data_log2 才是 CRDscore 要求的 log2 尺度表达矩阵
data_log2 <- as.data.frame(seurat_data / log(2))
exp = data_log2


target_genes <- c("HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "FDFT1", "SQLE",
                  "LSS", "CYP51A1", "TM7SF2", "MSMO1", "NSDHL", "HSD17B7", "DHCR24",
                  "EBP", "SC5D", "DHCR7")
target_genes = intersect(target_genes, rownames(pbmc1))

# --- 3. 计算评分 (方法 A: CRDscore) ---
score <- cal_CRDscore(expr = exp, n.bins = 50, circadians = 
                        target_genes, study.type = "scRNAseq")
gc()
score = as.data.frame(score)
score$id = rownames(score)

meta = pbmc1@meta.data
meta$id = rownames(meta) 
cli = read.csv("GSE178318_Cli.csv", header=T, check.names=F)

# 1. 按照患者 ID 合并临床信息
rt = merge(meta, cli, by.x= "orig.ident",by.y="Patient ID")

# 2. 按照细胞 ID 合并所有结果
rt1 = merge(score, rt, by.x = "id", by.y = "id")

# --- 绘图公共参数设置 ---
type_group = levels(factor(rt1$"type"))
my_comparisons = list()
if(length(type_group) >= 2){
  comp = combn(type_group, 2)
  for(i in 1:ncol(comp)){ my_comparisons[[i]] <- comp[,i] }
}

# ==================== 绘图 1: CRDscore ====================
data_p = rt1[,c('score','type')]
colnames(data_p) = c("score","Type")

p = ggplot(data_p, aes(x = Type, y = score, color = Type)) +
  stat_boxplot(geom="errorbar", width=0.6) +
  geom_boxplot(alpha=0.7, outlier.shape = NA, size=0.7, width=0.7, fatten=0.7) +
  theme_bw() + theme(panel.grid=element_blank()) +
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") +
  ggtitle(paste0(file_base_name, "+GSE178318+CRDscore")) +
  ylim(-0.2, 0.2) + ylab("CRDScore") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 15, face = "bold", color = "black"),
        axis.text.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.y = element_text(size = 15, face = "bold", color = "black"),
        axis.title.x = element_blank(), legend.position = "none",
        plot.title = element_text(face = "bold", size = 15, hjust = 0.5))

print(p)
ggsave(file.path(PLOTS_DIR, paste0(file_base_name, "_GSE178318_CRDscore.png")), plot = p, width = 7, height = 6, dpi = 300)
