# ==============================================================================
# Script 4: Monocytes 亚群精细分析 (Leiden 聚类 + 全局坐标投影 + 差异基因导出)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)
library(patchwork)
library(reticulate)

# 确保 Leiden 环境
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)

# 1. 读取上一步处理好的 PBMC 数据
if(file.exists('pbmc_corrected.qs')){
  pbmc <- qread('pbmc_corrected.qs')
} else {
  sc_by_tissue <- qread('sc_by_tissue_louvain.qs')
  pbmc <- sc_by_tissue[['PBMC']]
}

# ------------------------------------------------------------------------------
# 2. 提取 Monocytes 亚群并运行 Leiden 聚类
# ------------------------------------------------------------------------------
print("🚀 步骤1: 提取 Monocytes 亚群并运行子集精细聚类...")

# 提取单核细胞
mono_cells <- subset(pbmc, subset = celltype %in% grep("Monocytes", unique(pbmc$celltype), value = TRUE))

# 子集精细降维聚类
mono_cells <- FindVariableFeatures(mono_cells, selection.method = "vst", nfeatures = 2000)
mono_cells <- ScaleData(mono_cells)
mono_cells <- RunPCA(mono_cells, verbose = FALSE)

# 使用 Leiden 算法 (algorithm = 4)
mono_cells <- FindNeighbors(mono_cells, dims = 1:15)
mono_cells <- FindClusters(mono_cells, resolution = 1.0, algorithm = 4) 

# 生成亚群专属编号，防止与总图的 0,1,2 混淆
mono_cells$mono_sub_clusters <- paste0("Mono_", mono_cells$seurat_clusters)

# ------------------------------------------------------------------------------
# 3. 差异基因分析 (直接导出 Leiden 亚群 Marker)
# ------------------------------------------------------------------------------
print("🚀 步骤2: 计算单核细胞亚群 Marker 基因...")

Idents(mono_cells) <- "mono_sub_clusters"

# 计算所有阳性 Marker
all_markers <- FindAllMarkers(mono_cells, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              verbose = FALSE)

# 直接取 Top 20
top20_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

# 导出 Marker 列表
write.csv(all_markers, file.path("files", "Markers_All_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
write.csv(top20_markers, file.path("files", "Markers_Top20_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
print("  ✅ 单核细胞亚群 Marker 列表已导出。")

# ------------------------------------------------------------------------------
# 4. 全局坐标绘图 (2x2 网格排版)
# ------------------------------------------------------------------------------
print("🚀 步骤3: 正在生成基于总图坐标的 2x2 网格 UMAP...")

# 将亚群信息映射回 pbmc 总对象
pbmc$mono_viz_group <- "Others"
pbmc@meta.data[Cells(mono_cells), "mono_viz_group"] <- as.character(mono_cells$mono_sub_clusters)

# 获取各温度组下的单核细胞 Barcode
cells_cold <- Cells(subset(mono_cells, subset = Group == "Cold_4C"))
cells_rt   <- Cells(subset(mono_cells, subset = Group == "RT_25C"))
cells_tn   <- Cells(subset(mono_cells, subset = Group == "TN_30C"))

# ==========================================
# 绘图: 基于全局坐标的 Monocyte 亚群 (按温度分组)
# ==========================================
p_total_cls <- DimPlot(pbmc, cells = Cells(mono_cells), reduction = "umap", group.by = "mono_viz_group", label = TRUE, repel = TRUE) + 
  ggtitle("Monocyte Subgroups - Leiden (Total)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_cls <- DimPlot(pbmc, cells = cells_cold, reduction = "umap", group.by = "mono_viz_group", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_cls <- DimPlot(pbmc, cells = cells_rt, reduction = "umap", group.by = "mono_viz_group", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_cls <- DimPlot(pbmc, cells = cells_tn, reduction = "umap", group.by = "mono_viz_group", label = FALSE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 合并图层，应用 axes = "collect" 对齐边界
p_final_cls <- (p_total_cls | p_cold_cls) / (p_rt_cls | p_tn_cls) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

file_name_cls <- file.path("pictures", "UMAP_Grid_Monocytes_Subgroups_Leiden_Global.png")
ggsave(filename = file_name_cls, plot = p_final_cls, width = 15, height = 11, dpi = 300)

# 保存提取的单核细胞对象供后续独立调用
qsave(mono_cells, "pbmc_monocytes_sub-clustered.qs")

print("✅ 脚本 04 执行完毕，结果已保存！")