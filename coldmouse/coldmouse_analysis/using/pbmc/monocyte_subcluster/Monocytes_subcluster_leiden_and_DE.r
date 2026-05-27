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

# --- 【自定义注释与标签映射】 ---
# 定义 Seurat 默认聚类编号(0-5)对应的细胞注释
cluster_annotations <- c(
  "1" = "Activated/Mature Monocytes",
  "2" = "Pro-inflammatory classic Monocytes",
  "3" = "FKBP5+ classic Monocytes",
  "4" = "FKBP5+ non-classic Monocytes",
  "5" = "Antigen-Presenting Monocytes",
  "6" = "ISG high Monocytes"
)

# 1. 生成 Mono_1 到 Mono_6 的简洁编号作为主要聚类 ID
mono_cells$mono_cluster_id <- paste0("Mono_", as.numeric(as.character(mono_cells$seurat_clusters)))
cluster_number = sort(unique((mono_cells$seurat_clusters)))
# 固定 levels，确保顺序为 1 到 6
mono_cells$mono_cluster_id <- factor(mono_cells$mono_cluster_id, levels = paste0("Mono_", cluster_number))

# 2. 另开一列存储纯文本注释（备用）
mono_cells$mono_annotation <- unname(cluster_annotations[as.character(mono_cells$seurat_clusters)])

# 3. 提前定义好图例所需的 labels 映射向量 (供 scale_color_discrete 使用)
legend_labels <- paste0("Mono_", 1:6, ": ", cluster_annotations)
names(legend_labels) <- paste0("Mono_", 1:6)
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 3. 差异基因分析 (直接导出 Leiden 亚群 Marker)
# ------------------------------------------------------------------------------
print("🚀 步骤2: 计算单核细胞亚群 Marker 基因...")

Idents(mono_cells) <- "mono_cluster_id"

# 计算所有阳性 Marker
all_markers <- FindAllMarkers(mono_cells, 
                              only.pos = TRUE, 
                              min.pct = 0.25, 
                              logfc.threshold = 0.25,
                              verbose = FALSE)

# 直接取 Top 30
top20_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

# 导出 Marker 列表
write.csv(all_markers, file.path("files", "Markers_All_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
write.csv(top20_markers, file.path("files", "Markers_Top30_Leiden_Monocytes_Subgroups.csv"), row.names = FALSE)
print("  ✅ 单核细胞亚群 Marker 列表已导出。")
mono_cells = qread("pbmc_monocytes_sub-clustered.qs")
# ------------------------------------------------------------------------------
# 4. 局部坐标绘图 (2x2 网格排版，不映射回原图)
# ------------------------------------------------------------------------------
print("🚀 步骤3: 正在生成基于子集局部坐标的 2x2 网格 UMAP...")

# ==============================================================================
# 【新增步骤】自动获取全局 UMAP 坐标范围，并增加 5% 的边距（防止边缘细胞被切掉）
# ==============================================================================
umap_coords <- Embeddings(mono_cells, reduction = "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])

# 适当增加一点 buffer（边距）
x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05
global_xlim <- x_range + c(-x_buffer, x_buffer)
global_ylim <- y_range + c(-y_buffer, y_buffer)

# ==========================================
# 绘图 A: 基于局部坐标的 Monocyte 亚群 (直接画 mono_cells)
# ==========================================
p_total_local <- DimPlot(mono_cells, reduction = "umap", group.by = "mono_cluster_id", label = TRUE, repel = TRUE) + 
  scale_color_discrete(labels = legend_labels) + 
  ggtitle("Monocyte Subgroups - Leiden (Local Total)") + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_local <- DimPlot(subset(mono_cells, subset = Group == "Cold_4C"), reduction = "umap", group.by = "mono_cluster_id", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_local <- DimPlot(subset(mono_cells, subset = Group == "RT_25C"), reduction = "umap", group.by = "mono_cluster_id", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_local <- DimPlot(subset(mono_cells, subset = Group == "TN_30C"), reduction = "umap", group.by = "mono_cluster_id", label = FALSE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# ==========================================
# 合并图层并统一坐标轴
# ==========================================
# 使用 & coord_cartesian 可以把相同的坐标轴范围强加给拼图中的每一个子图
p_final_local <- (p_total_local | p_cold_local) / (p_rt_local | p_tn_local) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9)) & 
  coord_cartesian(xlim = global_xlim, ylim = global_ylim) # 【重点修改】在这里统一坐标轴范围

file_name_png <- file.path("pictures", "UMAP_Grid_Monocytes_Subgroups_Leiden_Local.png")
file_name_pdf <- file.path("pictures", "UMAP_Grid_Monocytes_Subgroups_Leiden_Local.pdf")
ggsave(filename = file_name_png, plot = p_final_local, width = 15, height = 11, dpi = 300)
ggsave(filename = file_name_pdf, plot = p_final_local, width = 15, height = 11)
# ------------------------------------------------------------------------------
# 5. 全局坐标绘图 (2x2 网格排版，映射回原图)
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在生成基于总图坐标的 2x2 网格 UMAP...")

# 将亚群信息映射回 pbmc 总对象
pbmc$mono_viz_group <- "Others"
pbmc@meta.data[Cells(mono_cells), "mono_viz_group"] <- as.character(mono_cells$mono_cluster_id)

# 固定总图中图例的顺序，并将 Others 放最后
pbmc$mono_viz_group <- factor(pbmc$mono_viz_group, levels = c(paste0("Mono_", 1:6), "Others"))

# 为全局图例追加 "Others" 的映射
legend_labels_global <- c(legend_labels, "Others" = "Others")

# 获取各温度组下的单核细胞 Barcode
cells_cold <- Cells(subset(mono_cells, subset = Group == "Cold_4C"))
cells_rt   <- Cells(subset(mono_cells, subset = Group == "RT_25C"))
cells_tn   <- Cells(subset(mono_cells, subset = Group == "TN_30C"))

# ==========================================
# 绘图 B: 基于全局坐标的 Monocyte 亚群 (按温度分组)
# ==========================================
p_total_cls <- DimPlot(pbmc, cells = Cells(mono_cells), reduction = "umap", group.by = "mono_viz_group", label = TRUE, repel = TRUE) + 
  scale_color_discrete(labels = legend_labels_global) + # 【重点修改】替换全局图例文字
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
# 保存提取的单核细胞对象供后续独立调用 (该对象现在自带局部的 UMAP 坐标和完整的注释信息)
qsave(mono_cells, "pbmc_monocytes_sub-clustered.qs")

print("✅ 脚本 04 执行完毕，局部/全局 2x2 UMAP 结果均已保存！")

# ------------------------------------------------------------------------------
# 6. Monocytes 组间差异基因分析 (Cold_4C vs RT_25C)
# ------------------------------------------------------------------------------
print("🚀 步骤5: 正在计算 Monocytes 在 Cold_4C 和 RT_25C 之间的差异表达基因 (DEG)...")

# 将细胞的身份标识 (Idents) 切换为实验组别 (Group)
Idents(mono_cells) <- "Group"

# 计算 Cold_4C 对比 RT_25C 的差异基因 (ident.1 vs ident.2)
deg_cold_vs_rt <- FindMarkers(mono_cells, 
                              ident.1 = "Cold_4C", 
                              ident.2 = "RT_25C",
                              logfc.threshold = 0.25,
                              min.pct = 0.1,  
                              verbose = FALSE)

# 提取行名(基因名)为单独一列，并按照 avg_log2FC 进行降序排列
deg_cold_vs_rt <- deg_cold_vs_rt %>% 
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC)) 

# 导出组间差异基因列表
deg_filename <- file.path("files", "DEG_Monocytes_Cold_4C_vs_RT_25C.csv")
write.csv(deg_cold_vs_rt, deg_filename, row.names = FALSE)

print(paste("  ✅ 组间差异基因 (Cold_4C vs RT_25C) 计算完成，已导出至:", deg_filename))
# ------------------------------------------------------------------------------
# 7. seurat cluster 3 / 4 分别与其他细胞做差异基因分析
# ------------------------------------------------------------------------------
print("🚀 步骤6: 正在计算 seurat cluster 3 和 4 分别 vs 其他细胞的差异基因...")

# 确保输出目录存在
dir.create("files", showWarnings = FALSE, recursive = TRUE)
mono_cells = qread('pbmc_monocytes_sub-clustered.qs')
# 将身份标识切换回 seurat_clusters
mono_cells$seurat_clusters <- factor(as.character(mono_cells$seurat_clusters))
Idents(mono_cells) <- "mono_cluster_id"

target_clusters <- c("Mono_3", "Mono_4")

# 检查 cluster 是否存在
missing_clusters <- setdiff(target_clusters, levels(Idents(mono_cells)))
if(length(missing_clusters) > 0){
  stop(paste0("以下 seurat cluster 不存在: ", paste(missing_clusters, collapse = ", ")))
}

deg_list <- list()
top20_list <- list()

for(cl in target_clusters){
  
  print(paste0("  🔍 正在计算 cluster ", cl, " vs Others ..."))
  
  # ident.2 = NULL 表示 ident.1 与所有其他细胞比较
  deg_tmp <- FindMarkers(
    mono_cells,
    ident.1 = cl,
    ident.2 = NULL,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  deg_tmp <- deg_tmp %>%
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC)) %>%
    mutate(
      cluster = paste0("cluster_", cl),
      comparison = paste0("cluster_", cl, "_vs_others")
    )
  
  # 按 avg_log2FC 取前 20 个上调差异基因
  top20_tmp <- deg_tmp %>%
    slice_max(n = 20, order_by = avg_log2FC, with_ties = FALSE)
  
  # 分别导出完整 DEG 和 Top20 DEG
  all_filename <- file.path("files", paste0("DEG_monoCluster_", cl, "_vs_Others_All.csv"))
  top20_filename <- file.path("files", paste0("DEG_monoCluster_", cl, "_vs_Others_Top20.csv"))
  
  write.csv(deg_tmp, all_filename, row.names = FALSE)
  write.csv(top20_tmp, top20_filename, row.names = FALSE)
  
  deg_list[[cl]] <- deg_tmp
  top20_list[[cl]] <- top20_tmp
  
  print(paste0("  ✅ seurat cluster ", cl, " vs Others 完成，Top20 已导出至: ", top20_filename))
}

# 合并 cluster 3 和 4 的 Top20 结果，方便查看
top20_combined <- bind_rows(top20_list)

combined_filename <- file.path("files", "DEG_monoCluster_3_4_vs_Others_Top20_Combined.csv")
write.csv(top20_combined, combined_filename, row.names = FALSE)

print(paste0("✅ mono cluster 3 和 4 vs Others 的 Top20 差异基因合并表已导出至: ", combined_filename))