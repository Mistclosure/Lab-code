# ==============================================================================
# 单细胞 WGCNA 分析流程：融合优化版参数
# 结果统一输出至 WGCNA_variable_genes 文件夹
# ==============================================================================

library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
library(WGCNA)
library(reshape2)
library(stringr)
library(qs) 

options(stringsAsFactors = FALSE) 

# 设置主工作目录（原始数据所在位置）
base_dir <- '/mnt/disk1/qiuzerui/downloads/CRC/GSE132465_GSE231559'
setwd(base_dir)
output_dir <- "WGCNA_variable_genes"
if (!dir.exists(output_dir)) dir.create(output_dir)

# ------------------------------------------------------------------------------
# 1. 读取数据与提取表达矩阵
# ------------------------------------------------------------------------------
cat("1. 正在读取单细胞数据...\n")
sc_obj <- qread("Merged_malig.qs")
sc_obj = JoinLayers(sc_obj)
# 确保从 data 图层提取标准化后的表达量
expr_mat <- as.matrix(GetAssayData(sc_obj, assay = "RNA", layer = "data"))

cell_meta <- data.frame(
  CellID  = colnames(sc_obj),
  Cluster = as.factor(sc_obj$seurat_clusters)
)
rownames(cell_meta) <- cell_meta$CellID

# 切换到输出文件夹
setwd(paste0(base_dir, "/", output_dir))
cat(paste("▶ 工作目录已切换至输出文件夹。\n"))

# ------------------------------------------------------------------------------
# 2. 伪细胞构建 (参考你的逻辑，按 cluster 划分)
# ------------------------------------------------------------------------------
cat("2. 正在构建伪细胞 (Pseudocells)...\n")
pseudocell.size <- 30
pseudo_df_list <- list() 

for (cluster_id in levels(cell_meta$Cluster)) {
  cluster_cells <- rownames(cell_meta[cell_meta$Cluster == cluster_id, ])
  if (length(cluster_cells) == 0) next
  
  pseudo_ids <- floor((seq_along(cluster_cells) - 1) / pseudocell.size)
  pseudo_ids <- paste0("Cluster", cluster_id, "_Pseudo", pseudo_ids)
  
  pseudo_df_list[[cluster_id]] <- data.frame(
    CellID = sample(cluster_cells), 
    PseudoCell = pseudo_ids,
    stringsAsFactors = FALSE
  )
}

pseudo_df <- do.call(rbind, pseudo_df_list)
expr_mat_t <- t(expr_mat[, pseudo_df$CellID])

cat("   正在聚合伪细胞平均表达量...\n")
pseudocell_expr_t <- aggregate(expr_mat_t, by = list(PseudoCell = pseudo_df$PseudoCell), FUN = mean)
rownames(pseudocell_expr_t) <- pseudocell_expr_t$PseudoCell
pseudocell_expr <- t(pseudocell_expr_t[, -1])

# ------------------------------------------------------------------------------
# 3. 提取高变基因 (Top 3000) 作为 WGCNA 输入
# ------------------------------------------------------------------------------
cat("3. 正在提取 Top 高变基因...\n")
expressed_ratio <- rowSums(pseudocell_expr > 0) / ncol(pseudocell_expr)
keep_genes_basic <- expressed_ratio > 0.05
pseudocell_expr_filtered <- pseudocell_expr[keep_genes_basic, ]

gene_var <- apply(pseudocell_expr_filtered, 1, var)
top_hvgs <- names(sort(gene_var, decreasing = TRUE))[1:3000]
target_expr_mat <- pseudocell_expr_filtered[top_hvgs, ]

# ------------------------------------------------------------------------------
# 4. WGCNA 网络构建 (融合参考代码参数)
# ------------------------------------------------------------------------------
cat("4. 正在进行 WGCNA 网络构建与参数计算...\n")
allowWGCNAThreads()

dataExpr <- as.matrix(t(target_expr_mat))

gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK) dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]

powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
# 参考代码使用了 signed 评估，但后续构建是 unsigned。这里统一使用 unsigned 评估更严谨
sft = pickSoftThreshold(dataExpr, powerVector = powers, networkType = "signed", verbose = 5)

png(filename = "WGCNA_SoftThreshold.png", width = 9, height = 5, units = "in", res = 300)
par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], labels = powers, col = "red")
abline(h = 0.85, col = "red")

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
dev.off()

# 容错处理：如果没有选出最优值，强制使用 12 或更高
power <- sft$powerEstimate
if(is.na(power) | power < 12) power <- 14 
cat(paste("使用的软阈值 (power) 为:", power, "\n"))

# --- 核心网络构建参数调整 ---
corType = "bicor"         # 强烈推荐保留 bicor 以压制单细胞极端值
maxPOutliers = 0.05       # 配合 bicor
cor <- WGCNA::cor         # 防止与其他包的 cor 函数冲突

net = blockwiseModules(dataExpr, power = power, maxBlockSize = ncol(dataExpr),
                       TOMType = "signed", 
                       minModuleSize = 10,       # 【参考代码调整】允许更小的精细模块
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,    
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE, corType = corType,
                       maxPOutliers = maxPOutliers, loadTOMs = TRUE,
                       saveTOMFileBase = "dataExpr.tom", verbose = 3)

# ------------------------------------------------------------------------------
# 5. 模块导出与 TOM 网络热图绘制
# ------------------------------------------------------------------------------
cat("5. 正在导出模块并绘制 TOM 热图...\n")
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)

# 导出基因与模块的对应关系 (参考代码格式)
Labels_out = data.frame(Gene = colnames(dataExpr), ModuleLabel = moduleLabels, ModuleColor = moduleColors)
write.table(Labels_out, file="WGCNA_moduleLabels.txt", sep="\t", quote=F, row.names=F)

# 绘制聚类树和模块颜色
png(filename = "WGCNA_Dendro_Colors.png", width = 8, height = 5, units = "in", res = 300)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# TOM网络热图绘制 (参考代码使用了 power=7)
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1 - TOM
plotTOM = dissTOM^15 # 【参考代码调整】提高幂次，让热图结构更清晰
diag(plotTOM) = NA

png(filename = "WGCNA_Network_Heatmap.png", width = 9, height = 9, units = "in", res = 300)
TOMplot(plotTOM, net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
        main = "Network heatmap plot",
        col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()

# ------------------------------------------------------------------------------
# 6. 临床特征关联 (Metastasis)
# ------------------------------------------------------------------------------
cat("6. 正在计算模块与临床特征的关联...\n")
MEs = net$MEs
colnames(MEs) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs = orderMEs(MEs)

# 提取 Metastasis 表型并与伪细胞对齐
raw_met <- sc_obj@meta.data[pseudo_df$CellID, "metastasis"]
cell_metastasis <- ifelse(is.character(raw_met) | is.factor(raw_met), as.numeric(as.factor(raw_met)) - 1, as.numeric(raw_met))

cell_traits_df <- data.frame(PseudoCell = pseudo_df$PseudoCell, metastasis = cell_metastasis)
pseudo_traits <- aggregate(metastasis ~ PseudoCell, data = cell_traits_df, FUN = mean)
rownames(pseudo_traits) <- pseudo_traits$PseudoCell

common_samples <- intersect(rownames(MEs), rownames(pseudo_traits))
MEs <- MEs[common_samples, ]
cli <- pseudo_traits[common_samples, "metastasis", drop = FALSE]

moduleTraitCor <- cor(MEs, cli, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, length(common_samples))

png(filename = "WGCNA_Module_Trait_Heatmap.png", width = 5, height = 7, units = "in", res = 300)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
par(mar = c(10, 6, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Metastasis",
               yLabels = names(MEs), ySymbols = names(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMatrix, setStdMargins = FALSE,
               cex.text = 0.65, zlim = c(-1, 1),
               main = "Module-trait relationships")
dev.off()

cat("\n✅ 分析已完成！请查看生成的 WGCNA_Network_Heatmap.png 以确认模块区分度。\n")