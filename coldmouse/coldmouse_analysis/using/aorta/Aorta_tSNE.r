# ==============================================================================
# Script 02.01: Aorta 独立降维聚类 (tSNE版) 与全局 SingleR 注释
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
# 1. 读取数据并提取 Aorta 组织
# ------------------------------------------------------------------------------
print("🚀 步骤1: 读取 sc_by_tissue 数据并提取 Aorta...")
# 读取脚本2中第5步或第6步保存下来的 qs 文件
sc_by_tissue <- qread(paste0('sc_by_tissue_', cluster_method, '.qs'))

if (!"Aorta" %in% names(sc_by_tissue)) {
  stop("❌ 错误：sc_by_tissue 数据中未找到 Aorta 组织，请检查输入文件！")
}

aorta_obj <- sc_by_tissue[["Aorta"]]
print(paste(">>> 成功提取 Aorta，当前细胞数:", ncol(aorta_obj)))

# ------------------------------------------------------------------------------
# 4. 重新进行降维聚类 (使用 tSNE 替代 UMAP)
# ------------------------------------------------------------------------------
print(paste0("🚀 步骤4: 重新运行标准化与降维聚类 (使用 tSNE)..."))

# 确保 Seurat v5 layer 合并
aorta_obj <- JoinLayers(aorta_obj)

# 重新走一遍标准流程
aorta_obj <- NormalizeData(aorta_obj, normalization.method = "LogNormalize", scale.factor = 10000)
aorta_obj <- FindVariableFeatures(aorta_obj, selection.method = "vst")
aorta_obj <- ScaleData(aorta_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
aorta_obj <- RunPCA(aorta_obj, npcs = 50, verbose = FALSE)

# 聚类
aorta_obj <- FindNeighbors(aorta_obj, reduction = "pca", dims = 1:50)
aorta_obj <- FindClusters(aorta_obj, resolution = 1.0, algorithm = algo_id) 

# 【关键修改】：使用 tSNE 进行降维
print("   🔎 正在运行 tSNE 降维...")
aorta_obj <- RunTSNE(aorta_obj, reduction = "pca", dims = 1:50)

# ------------------------------------------------------------------------------
# 5. 重新运行 SingleR 注释
# ------------------------------------------------------------------------------
print("🚀 步骤5: 运行 SingleR 对 Aorta 细胞进行全局标注...")
mouse_ref <- celldex::MouseRNAseqData()

expr_mat <- GetAssayData(aorta_obj, assay = "RNA", layer = "data")
pred_res <- SingleR(test = expr_mat, ref = mouse_ref, labels = mouse_ref$label.main)
aorta_obj$SingleR.labels <- pred_res$labels
print("   ✅ Aorta 全局 SingleR 注释完成")

# ------------------------------------------------------------------------------
# 6. 绘制细胞聚类图片 (分温度) 
# ------------------------------------------------------------------------------
print("🚀 步骤6: 绘制 Aorta 的 tSNE 聚类图并输出...")

# 过滤极小亚群
min_cells_threshold <- 20
cell_counts <- table(aorta_obj$SingleR.labels)
keep_labels <- names(cell_counts[cell_counts >= min_cells_threshold])

if (length(keep_labels) == 0) {
  stop("   ⚠️ Aorta 中没有符合条件的亚群。")
}

aorta_obj <- subset(aorta_obj, subset = SingleR.labels %in% keep_labels)

# 规范化 label 因子水平
all_labels <- sort(unique(aorta_obj$SingleR.labels))
aorta_obj$SingleR.labels <- factor(aorta_obj$SingleR.labels, levels = all_labels)

# ==========================================
# 绘图: 基于 SingleR 注释 (按 Annotation 分组)
# 注意：这里已将 reduction 全部改为 "tsne"
# ==========================================
p_total_anno <- DimPlot(aorta_obj, reduction = "tsne", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
  ggtitle(paste("Aorta -", cluster_method, "(Annotation - tSNE)")) + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_anno <- DimPlot(subset(aorta_obj, subset = Group == "Cold_4C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_anno <- DimPlot(subset(aorta_obj, subset = Group == "RT_25C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_anno <- DimPlot(subset(aorta_obj, subset = Group == "TN_30C"), reduction = "tsne", group.by = "SingleR.labels", label = FALSE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 拼接成 2x2 网格
p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

# 导出图片
file_name_anno <- file.path("pictures", paste0("tSNE_Grid_Annotation_", cluster_method, "_Aorta.png"))
ggsave(filename = file_name_anno, plot = p_final_anno, width = 15, height = 11, dpi = 300)

print(paste("   ✅ Aorta 的 Annotation 版 tSNE 图像已保存至:", file_name_anno))

# (可选) 独立保存一下包含 tSNE 的 Aorta 专属对象，方便后续探索
#qsave(aorta_obj, file.path("files", paste0('sc_Aorta_tSNE_', cluster_method, '.qs')))
#print("   ✅ 独立处理的 Aorta 数据已保存为 qs 文件。")