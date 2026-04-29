# ==============================================================================
# Script 02.02: PBMC 独立降维聚类 (tSNE版) 与全局 SingleR 注释
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 确保输出文件夹存在
dir.create("pictures", showWarnings = FALSE)
dir.create("files", showWarnings = FALSE)

# 【参数配置】保持与之前一致的聚类方法
cluster_method <- "louvain" 
algo_id <- ifelse(cluster_method == "leiden", 4, 1)

library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(SingleR)
library(celldex)

# ------------------------------------------------------------------------------
# 1. 读取 PBMC 数据
# ------------------------------------------------------------------------------
print("🚀 步骤1: 读取 pbmc_corrected.qs 数据...")
# 直接读取仅包含 PBMC 的对象
pbmc_obj <- qread("pbmc_corrected.qs")

print(paste(">>> 成功读取 PBMC 数据，当前细胞数:", ncol(pbmc_obj)))

# ------------------------------------------------------------------------------
# 4. 重新进行降维聚类 (使用 tSNE 替代 UMAP)
# ------------------------------------------------------------------------------
print(paste0("🚀 步骤4: 重新运行标准化与降维聚类 (使用 tSNE)..."))

# 确保 Seurat v5 layer 合并
pbmc_obj <- JoinLayers(pbmc_obj)

# 重新走一遍标准流程
pbmc_obj <- NormalizeData(pbmc_obj, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc_obj <- FindVariableFeatures(pbmc_obj, selection.method = "vst")
pbmc_obj <- ScaleData(pbmc_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
pbmc_obj <- RunPCA(pbmc_obj, npcs = 50, verbose = FALSE)

# 聚类
pbmc_obj <- FindNeighbors(pbmc_obj, reduction = "pca", dims = 1:50)
pbmc_obj <- FindClusters(pbmc_obj, resolution = 1.0, algorithm = algo_id) 

# 【关键执行】：使用 tSNE 进行降维
print("   🔎 正在运行 tSNE 降维...")
pbmc_obj <- RunTSNE(pbmc_obj, reduction = "pca", dims = 1:50)

# ------------------------------------------------------------------------------
# 5. 重新运行 SingleR 注释
# ------------------------------------------------------------------------------
print("🚀 步骤5: 运行 SingleR 对 PBMC 细胞进行全局标注...")
mouse_ref <- celldex::MouseRNAseqData()

expr_mat <- GetAssayData(pbmc_obj, assay = "RNA", layer = "data")
pred_res <- SingleR(test = expr_mat, ref = mouse_ref, labels = mouse_ref$label.main)
pbmc_obj$SingleR.labels <- pred_res$labels
print("   ✅ PBMC 全局 SingleR 注释完成")

# ------------------------------------------------------------------------------
# 6. 绘制细胞聚类图片 (分温度) 
# ------------------------------------------------------------------------------
print("🚀 步骤6: 绘制 PBMC 的 tSNE 聚类图并输出...")

# 过滤极小亚群
min_cells_threshold <- 20
cell_counts <- table(pbmc_obj$SingleR.labels)
keep_labels <- names(cell_counts[cell_counts >= min_cells_threshold])

if (length(keep_labels) == 0) {
  stop("   ⚠️ PBMC 中没有符合条件的亚群。")
}

pbmc_obj <- subset(pbmc_obj, subset = SingleR.labels %in% keep_labels)

# 规范化 label 因子水平
all_labels <- sort(unique(pbmc_obj$SingleR.labels))
pbmc_obj$SingleR.labels <- factor(pbmc_obj$SingleR.labels, levels = all_labels)

# ==========================================
# 绘图: 基于 SingleR 注释 (按 Annotation 分组)
# 注意：基于 tSNE 的绘图逻辑
# ==========================================
p_total_anno <- DimPlot(pbmc_obj, reduction = "tsne", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle(paste("PBMC -", cluster_method, "(Annotation - tSNE)")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_anno <- DimPlot(subset(pbmc_obj, subset = Group == "Cold_4C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_anno <- DimPlot(subset(pbmc_obj, subset = Group == "RT_25C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_anno <- DimPlot(subset(pbmc_obj, subset = Group == "TN_30C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 拼接成 2x2 网格
p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

# 导出图片
file_name_anno <- file.path("pictures", paste0("tSNE_Grid_Annotation_", cluster_method, "_PBMC.png"))
ggsave(filename = file_name_anno, plot = p_final_anno, width = 15, height = 11, dpi = 300)

print(paste("   ✅ PBMC 的 Annotation 版 tSNE 图像已保存至:", file_name_anno))

# (可选) 如果你后续还需要用到这个含有 tSNE 和 SingleR 结果的 PBMC 对象，可以取消注释保存它
# qsave(pbmc_obj, file.path("files", paste0('sc_PBMC_tSNE_', cluster_method, '.qs')))
# print("   ✅ 独立处理的 PBMC 数据已保存为 qs 文件。")