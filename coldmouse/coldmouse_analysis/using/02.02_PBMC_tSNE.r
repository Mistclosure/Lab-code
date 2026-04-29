# ==============================================================================
# Script 02.02: PBMC 独立降维可视化 (极致精简版 - 仅运行 tSNE 并绘图)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 确保输出文件夹存在
dir.create("pictures", showWarnings = FALSE)
dir.create("files", showWarnings = FALSE)

# 【参数配置】
cluster_method <- "louvain" 

library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)

# ------------------------------------------------------------------------------
# 1. 读取 PBMC 数据 (已包含 PCA 和 SingleR 注释)
# ------------------------------------------------------------------------------
print("🚀 步骤1: 读取已处理好的 pbmc_obj 数据...")
sc_by_tissue = qread('sc_by_tissue_louvain.qs')
pbmc_obj = sc_by_tissue[['PBMC']]

print(paste(">>> 成功读取 PBMC 数据，当前细胞数:", ncol(pbmc_obj)))

# ------------------------------------------------------------------------------
# 2. 仅运行 tSNE 降维
# ------------------------------------------------------------------------------
print("🚀 步骤2: 基于已有的 PCA 结果计算 tSNE 坐标...")
# 直接调用已有的前 50 个主成分进行 tSNE 降维
pbmc_obj <- RunTSNE(pbmc_obj, reduction = "pca", dims = 1:50)

# ------------------------------------------------------------------------------
# 3. 绘制细胞 tSNE 图片 (分温度) 
# ------------------------------------------------------------------------------
print("🚀 步骤3: 绘制 PBMC 的 tSNE 分布图并输出...")

# 基于已有的 SingleR.labels 过滤极小亚群
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
# 绘图: 基于已有的 SingleR 注释 (按 Annotation 分组)
# 注意：使用 reduction = "tsne"
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

# (可选) 导出带有 tSNE 坐标的对象，以备后续其他分析需要
# qsave(pbmc_obj, file.path("files", paste0('sc_PBMC_tSNE_', cluster_method, '.qs')))