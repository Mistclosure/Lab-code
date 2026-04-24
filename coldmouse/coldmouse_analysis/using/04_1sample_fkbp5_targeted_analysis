# ==============================================================================
# Script 4: Monocytes 亚群精细分析 (Leiden 聚类 + 全局与局部坐标投影 + 差异基因导出)
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
mono_cells <- FindClusters(mono_cells, resolution = 0.35, algorithm = 4) 

# 【新增】运行局部的 UMAP 降维，为后续的局部坐标绘图提供数据
mono_cells <- RunUMAP(mono_cells, dims = 1:15, verbose = FALSE)

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
  slice_max(n = 30, order_by = avg_log2FC)

# 导出 Marker 列表
write.csv(all_markers, file.path("files", "Markers_All_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
write.csv(top20_markers, file.path("files", "Markers_Top30_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
print("  ✅ 单核细胞亚群 Marker 列表已导出。")

# ------------------------------------------------------------------------------
# 4. 局部坐标绘图 (2x2 网格排版，不映射回原图)
# ------------------------------------------------------------------------------
print("🚀 步骤3: 正在生成基于子集局部坐标的 2x2 网格 UMAP...")

# ==========================================
# 绘图 A: 基于局部坐标的 Monocyte 亚群 (直接画 mono_cells)
# ==========================================
p_total_local <- DimPlot(mono_cells, reduction = "umap", group.by = "mono_sub_clusters", label = TRUE, repel = TRUE) + 
  ggtitle("Monocyte Subgroups - Leiden (Local Total)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_local <- DimPlot(subset(mono_cells, subset = Group == "Cold_4C"), reduction = "umap", group.by = "mono_sub_clusters", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_local <- DimPlot(subset(mono_cells, subset = Group == "RT_25C"), reduction = "umap", group.by = "mono_sub_clusters", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_local <- DimPlot(subset(mono_cells, subset = Group == "TN_30C"), reduction = "umap", group.by = "mono_sub_clusters", label = FALSE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 合并图层，应用 axes = "collect" 对齐边界
p_final_local <- (p_total_local | p_cold_local) / (p_rt_local | p_tn_local) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

file_name_local <- file.path("pictures", "UMAP_Grid_Monocytes_Subgroups_Leiden_Local.png")
ggsave(filename = file_name_local, plot = p_final_local, width = 15, height = 11, dpi = 300)

# ------------------------------------------------------------------------------
# 5. 全局坐标绘图 (2x2 网格排版，映射回原图)
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在生成基于总图坐标的 2x2 网格 UMAP...")

# 将亚群信息映射回 pbmc 总对象
pbmc$mono_viz_group <- "Others"
pbmc@meta.data[Cells(mono_cells), "mono_viz_group"] <- as.character(mono_cells$mono_sub_clusters)

# 获取各温度组下的单核细胞 Barcode
cells_cold <- Cells(subset(mono_cells, subset = Group == "Cold_4C"))
cells_rt   <- Cells(subset(mono_cells, subset = Group == "RT_25C"))
cells_tn   <- Cells(subset(mono_cells, subset = Group == "TN_30C"))

# ==========================================
# 绘图 B: 基于全局坐标的 Monocyte 亚群 (按温度分组)
# ==========================================
p_total_cls <- DimPlot(pbmc, cells = Cells(mono_cells), reduction = "umap", group.by = "mono_viz_group", label = TRUE, repel = TRUE) + 
  ggtitle("Monocyte Subgroups - Leiden (Global Total)") + 
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

# 保存提取的单核细胞对象供后续独立调用 (该对象现在自带局部的 UMAP 坐标)
qsave(mono_cells, "pbmc_monocytes_sub-clustered.qs")

print("✅ 脚本 04 执行完毕，局部/全局 2x2 UMAP 结果均已保存！")

# ------------------------------------------------------------------------------
# 6. Monocytes 组间差异基因分析 (Cold_4C vs RT_25C)
# ------------------------------------------------------------------------------
print("🚀 步骤5: 正在计算 Monocytes 在 Cold_4C 和 RT_25C 之间的差异表达基因 (DEG)...")

# 将细胞的身份标识 (Idents) 切换为实验组别 (Group)
Idents(mono_cells) <- "Group"

# 计算 Cold_4C 对比 RT_25C 的差异基因 (ident.1 vs ident.2)
# 正的 log2FC 表示在 Cold_4C 中上调，负的表示在 Cold_4C 中下调
deg_cold_vs_rt <- FindMarkers(mono_cells, 
                              ident.1 = "Cold_4C", 
                              ident.2 = "RT_25C",
                              logfc.threshold = 0.25,
                              min.pct = 0.1,  # 稍微放宽 min.pct 以捕捉更多组间变化
                              verbose = FALSE)

# 提取行名(基因名)为单独一列，并按照 avg_log2FC 进行降序排列
deg_cold_vs_rt <- deg_cold_vs_rt %>% 
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) 

# 导出组间差异基因列表
deg_filename <- file.path("files", "DEG_Monocytes_Cold_4C_vs_RT_25C.csv")
write.csv(deg_cold_vs_rt, deg_filename, row.names = FALSE)

print(paste("  ✅ 组间差异基因 (Cold_4C vs RT_25C) 计算完成，已导出至:", deg_filename))

# ==============================================================================
# 【新增逻辑】: 提取 Fkbp5 表达状态与 2x2 网格 FeaturePlot
# ==============================================================================

# ------------------------------------------------------------------------------
# 7. 提取 Fkbp5 表达状态
# ------------------------------------------------------------------------------
print("🚀 步骤6: 正在提取 Fkbp5 表达状态...")

target_gene <- "Fkbp5"
tissue_name <- "PBMC"

fkbp5_counts <- GetAssayData(mono_cells, assay = "RNA", layer = "counts")[target_gene, ]
mono_cells$Fkbp5_Status <- ifelse(fkbp5_counts > 0, "Fkbp5+ Monocyte", "Fkbp5- Monocyte")
mono_cells$Fkbp5_Status <- factor(mono_cells$Fkbp5_Status, levels = c("Fkbp5- Monocyte", "Fkbp5+ Monocyte"))

print("  ✅ Fkbp5 表达状态提取完成。")

# ------------------------------------------------------------------------------
# 8. 生成 Fkbp5 表达量的 2x2 网格 FeaturePlot
# ------------------------------------------------------------------------------
print("🚀 步骤7: 正在绘制并保存 Fkbp5 表达量的 2x2 FeaturePlot...")

# 1. 提取全局 UMAP 的坐标极值，固定 X 和 Y 轴 (防止子图放大/变形)
umap_coords <- Embeddings(mono_cells, reduction = "umap")
x_lims <- range(umap_coords[, 1])
y_lims <- range(umap_coords[, 2])

# 2. 提取靶基因的最大表达量，固定颜色标尺 (防止低表达组出现“假红”)
max_expr <- max(FetchData(mono_cells, vars = target_gene)[, 1])

p_feat_total <- FeaturePlot(mono_cells, features = target_gene, reduction = "umap") + 
  scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
  xlim(x_lims) + ylim(y_lims) +
  ggtitle(paste(tissue_name, "Monocytes - Total")) +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_feat_cold <- FeaturePlot(subset(mono_cells, subset = Group == "Cold_4C"), features = target_gene, reduction = "umap") + 
  scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
  xlim(x_lims) + ylim(y_lims) +
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_feat_rt <- FeaturePlot(subset(mono_cells, subset = Group == "RT_25C"), features = target_gene, reduction = "umap") + 
  scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
  xlim(x_lims) + ylim(y_lims) +
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_feat_tn <- FeaturePlot(subset(mono_cells, subset = Group == "TN_30C"), features = target_gene, reduction = "umap") + 
  scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
  xlim(x_lims) + ylim(y_lims) +
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_feat_final <- (p_feat_total | p_feat_cold) / (p_feat_rt | p_feat_tn) + 
  plot_layout(guides = "collect", axes = "collect") 

# 保存拼接后的 FeaturePlot
ggsave(filename = paste0("pictures/Monocytes_Fkbp5_FeaturePlot_Grid_", tissue_name, ".png"), plot = p_feat_final, width = 15, height = 11, dpi = 300)

print(paste("  ✅ 2x2 FeaturePlot 已保存至: files/Monocytes_Fkbp5_FeaturePlot_Grid_", tissue_name, ".png"))
# ------------------------------------------------------------------------------
# 9. 生成 Fkbp5 状态的 2x2 网格 DimPlot (映射回全局坐标)
# ------------------------------------------------------------------------------
print("🚀 步骤8: 正在将 Fkbp5 状态映射回总图并生成全局 2x2 视图...")

# 将状态同步至全局对象
pbmc$Fkbp5_viz_status <- "Others"
pbmc@meta.data[Cells(mono_cells), "Fkbp5_viz_status"] <- as.character(mono_cells$Fkbp5_Status)

# 定义固定颜色
status_cols <- c("Fkbp5+ Monocyte" = "red", "Fkbp5- Monocyte" = "lightgrey")

# 基于总图坐标绘制
p_stat_total <- DimPlot(pbmc, cells = Cells(mono_cells), reduction = "umap", group.by = "Fkbp5_viz_status", cols = status_cols) + 
  ggtitle("Fkbp5 Status - Global Total") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_stat_cold <- DimPlot(pbmc, cells = cells_cold, reduction = "umap", group.by = "Fkbp5_viz_status", cols = status_cols) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_stat_rt <- DimPlot(pbmc, cells = cells_rt, reduction = "umap", group.by = "Fkbp5_viz_status", cols = status_cols) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_stat_tn <- DimPlot(pbmc, cells = cells_tn, reduction = "umap", group.by = "Fkbp5_viz_status", cols = status_cols) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_stat_final <- (p_stat_total | p_stat_cold) / (p_stat_rt | p_stat_tn) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 10))

file_name_stat_global <- file.path("pictures", paste0("UMAP_Grid_Monocytes_Fkbp5_Status_Global_", tissue_name, ".png"))
ggsave(filename = file_name_stat_global, plot = p_stat_final, width = 15, height = 11, dpi = 300)

print(paste("  ✅ 全局 Fkbp5 状态图已导出至:", file_name_stat_global))