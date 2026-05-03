# ==============================================================================
# 单细胞 WGCNA 分析流程 (恶性细胞子集伪细胞构建与共表达网络)
# ==============================================================================

# 1. 载入所需的 R 包 -----------------------------------------------------------
library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(SingleR)
library(WGCNA)
library(reshape2)
library(stringr)
library(qs)

setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465')

# 2. 读取数据与提取表达矩阵 ----------------------------------------------------
cat("正在加载恶性细胞 Seurat 对象...\n")
malig_seurat <- qread("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/Malignant_RNA_assay.qs")

# 提取 RNA assay 的 data 图层作为表达矩阵
malig_seurat = JoinLayers(malig_seurat)
expr_mat <- as.matrix(GetAssayData(malig_seurat, assay = "RNA", layer = "data"))

# 3. 提取与整理细胞聚类信息 ----------------------------------------------------
cat("正在整理细胞聚类信息...\n")
meta_data <- malig_seurat@meta.data

# 构建初步的细胞聚类数据框
cell_types_df <- cbind(rownames(meta_data), as.character(meta_data$seurat_clusters))
rownames(cell_types_df) <- rownames(meta_data)
colnames(cell_types_df) <- c("CellID", "Celltype")
cell_types_df <- as.data.frame(cell_types_df)

# 与原始 metadata 合并，确保信息对齐
set.seed(9)
merged_metadata <- merge(malig_seurat@meta.data, cell_types_df, by = 0)

# 提取核心的 CellID 与 Celltype 信息
cell_types_df <- merged_metadata[, c(1, ncol(merged_metadata))] 
colnames(cell_types_df) <- c("CellID", "Celltype")
rownames(cell_types_df) <- cell_types_df[, 1]
cell_types_df$Celltype <- as.factor(cell_types_df$Celltype)

# 根据整理好的 CellID 对表达矩阵进行切片对齐
expr_subset <- expr_mat[, cell_types_df$CellID]
expr_subset_mat <- as.matrix(expr_subset)

# ==============================================================================
# 4. 伪细胞 (Pseudo-cell) 构建 -------------------------------------------------
# ==============================================================================
cat("正在构建伪细胞矩阵...\n")
pseudo_cell_size <- 30
pseudo_ids_list <- list()
cluster_levels <- levels(cell_types_df$Celltype)

# 4.1 为每个 Cluster 分配伪细胞 ID
for (i in 1:length(cluster_levels)) {
  cluster_id <- cluster_levels[i]
  cluster_cells <- rownames(cell_types_df[cell_types_df$Celltype == cluster_id, ])
  
  # 按指定大小分配批次 ID
  pseudo_ids <- floor(seq_along(cluster_cells) / pseudo_cell_size)
  pseudo_ids <- paste0(cluster_id, "_Cell", pseudo_ids)
  names(pseudo_ids) <- sample(cluster_cells)
  
  pseudo_ids_list[[i]] <- pseudo_ids
}

# 4.2 整合所有伪细胞 ID
all_pseudo_ids <- unlist(pseudo_ids_list)
all_pseudo_ids <- as.data.frame(all_pseudo_ids)
pseudo_cell_counts <- table(all_pseudo_ids)
valid_cell_barcodes <- rownames(all_pseudo_ids)

# 4.3 聚合表达量 (计算每个伪细胞的均值)
expr_for_aggregation <- expr_mat[, as.character(valid_cell_barcodes)]
expr_for_aggregation <- t(expr_for_aggregation)

aggregated_pseudo_expr <- aggregate(expr_for_aggregation, list(name = all_pseudo_ids[, 1]), FUN = mean)
rownames(aggregated_pseudo_expr) <- aggregated_pseudo_expr$name
aggregated_pseudo_expr <- aggregated_pseudo_expr[, -1]

# 4.4 转置回常规表达矩阵格式 (基因 x 细胞)
pseudo_cell_counts_mat <- as.matrix(pseudo_cell_counts)
pseudo_expr_matrix <- t(aggregated_pseudo_expr)[, rownames(pseudo_cell_counts_mat)]
qsave(pseudo_expr_matrix,'pseudo_expr_matrix.qs')
# ==============================================================================
# 5. 目标基因过滤与 WGCNA 数据准备 ---------------------------------------------
# ==============================================================================
cat("正在读取目标基因列表并过滤表达矩阵...\n")
ciliahub_genes <- read.csv("/mnt/disk1/qiuzerui/downloads/signature/ciliahub_genes_list.csv", header = TRUE, check.names = FALSE, row.names = 1)

# 提取交集基因的伪细胞表达矩阵
filtered_expr_matrix <- pseudo_expr_matrix[intersect(rownames(ciliahub_genes), rownames(pseudo_expr_matrix)), ]

# WGCNA 参数设置 (转置为 样本 x 基因 格式)
type <- "unsigned"
corType <- "pearson"
corFnc <- ifelse(corType == "pearson", cor, bicor)
maxPOutliers <- ifelse(corType == "pearson", 1, 0.05)
robustY <- ifelse(corType == "pearson", TRUE, FALSE)
dataExpr <- as.matrix(t(filtered_expr_matrix))

# 检查缺失值与低质量样本/基因
gsg <- goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) { 
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")))
  dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]
}

nGenes <- ncol(dataExpr)
nSamples <- nrow(dataExpr)

# ==============================================================================
# 6. WGCNA 网络构建 ------------------------------------------------------------
# ==============================================================================
cat("正在计算软阈值...\n")
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, networkType = "signed", verbose = 5)

# 6.1 保存软阈值图
png(filename = "WGCNA_SoftThreshold.png", width = 10, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n", main = "Scale independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
dev.off()

# 获取推荐软阈值，并设置容错
power <- sft$powerEstimate
if(is.na(power)) power <- 12 
softPower <- power

# 6.2 一步法构建共表达网络
cat("正在构建 Blockwise Modules...\n")
cor <- WGCNA::cor
net <- blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                        TOMType = "unsigned", minModuleSize = 10,
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE, corType = corType,
                        maxPOutliers = maxPOutliers, loadTOMs = TRUE,
                        saveTOMFileBase = paste0("dataExpr", ".tom"),
                        verbose = 3)

# 模块颜色转换
moduleLabels <- net$colors
moduleColors <- labels2colors(moduleLabels)

# 6.3 保存模块聚类树
png(filename = "WGCNA_ModuleColors_Dendro.png", width = 8, height = 6, units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# 6.4 提取模块特征基因 (MEs) 并保存热图 1
MEs <- net$MEs
MEs_col <- MEs
colnames(MEs_col) <- paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs), "ME", ""))))
MEs_col <- orderMEs(MEs_col)

png(filename = "WGCNA_Eigengene_Adjacency1.png", width = 7, height = 7, units = "in", res = 300)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marDendro = c(3, 3, 2, 4),
                      marHeatmap = c(3, 4, 2, 2), plotDendrograms = TRUE, xLabelsAngle = 90)
dev.off()

# ==============================================================================
# 7. 模块可视化与 TOM 网络分析 -------------------------------------------------
# ==============================================================================
cat("正在进行模块表达与 TOM 网络可视化...\n")

MEList <- moduleEigengenes(dataExpr, colors = moduleColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

# 保存特征基因热图 2
png(filename = "WGCNA_Eigengene_Adjacency2.png", width = 7, height = 7, units = "in", res = 300)
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marHeatmap = c(3, 4, 2, 2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

# 7.1 针对 Blue 模块的可视化
which.module <- "blue"
ME <- MEs[, paste0("ME", which.module)]

png(filename = "WGCNA_Blue_Module_Plot.png", width = 8, height = 8, units = "in", res = 300)
par(mfrow = c(2, 1), mar = c(0, 4.1, 4, 2.05))
plotMat(t(scale(dataExpr[, moduleColors == which.module])),
        nrgcols = 30, rlabels = FALSE, rcols = which.module,
        main = which.module, cex.main = 2)
par(mar = c(2, 2.3, 0.5, 0.8))
barplot(ME, col = which.module, main = "", cex.main = 2,
        ylab = "eigengene expression", xlab = "array sample")
dev.off()

# 7.2 TOM 网络热图绘制
load(net$TOMFiles, verbose = TRUE)
TOM <- as.matrix(TOM)
dissTOM <- 1 - TOM
plotTOM <- dissTOM^7
diag(plotTOM) <- NA

png(filename = "WGCNA_Network_Heatmap.png", width = 9, height = 9, units = "in", res = 300)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot, all genes",
        col = gplots::colorpanel(250, 'red', "orange", 'lemonchiffon'))
dev.off()

# 保存模块标签
Labels <- as.data.frame(moduleLabels)
Colors <- as.data.frame(moduleColors)
Labels <- cbind(Labels, Colors)
Colors <- cbind(id = rownames(Labels), Labels)
write.table(Colors, file = "moduleLabels.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# ==============================================================================
# 8. 模块与临床特征关联分析 ----------------------------------------------------
# ==============================================================================
cat("正在分析模块与临床特征的关联...\n")
cli <- read.csv("Cli.csv", header = TRUE, check.names = FALSE)
rownames(cli) <- cli$ID

# 确保 MEs 与 cli 样本顺序一致
common_samples <- intersect(rownames(MEs), rownames(cli))
MEs <- MEs[common_samples, ]
cli <- cli[common_samples, ]

moduleTraitCor <- cor(MEs, cli, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, length(common_samples))

# 保存临床关联热图
png(filename = "Module-trait.png", width = 7, height = 7, units = "in", res = 300)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
par(mar = c(10, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(cli),
               yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.65, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

cat("全部流程运行完毕！\n")