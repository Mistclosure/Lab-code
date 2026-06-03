# ==============================================================================
# 06_pbmc_monocyte_subcluster_introns.R
# coldmouse_introns：PBMC Monocyte 亚群精细聚类 + DEG 分析
# ==============================================================================
# 输入：coldmouse_introns_sc_by_tissue_no_harmony_louvain.qs (PBMC 部分)
# 输出：monocyte 子集对象、UMAP、markers、Cold_vs_RT/TN DEG
# 参考：Monocytes_subcluster_leiden_and_DE.r、PBMC_Monocyte_DE.r
# ==============================================================================

PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))

config <- load_config()
setup_performance(config)

library(Seurat)
library(tidyverse)
library(qs)
library(patchwork)

# 输出目录
output_dir <- config$paths$output_dir
files_dir <- file.path(output_dir, "files")
plots_dir <- file.path(output_dir, "plots")
objects_dir <- file.path(output_dir, "results", "objects")
subcluster_method <- "leiden"
subcluster_algorithm <- 4
subcluster_resolution <- 0.35
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)

# 读取按组织对象
sc_by_tissue_file <- file.path(output_dir,
                               paste0("coldmouse_introns_sc_by_tissue_no_harmony_",
                                      config$clustering$method, ".qs"))
sc_by_tissue <- qread(sc_by_tissue_file)
pbmc <- sc_by_tissue[["PBMC"]]

cat("=== PBMC Monocyte 亚群分析 ===\n")
cat("PBMC 总细胞数:", ncol(pbmc), "\n")
cat("SingleR 标签:", paste(sort(unique(pbmc$SingleR.labels)), collapse=", "), "\n")

# ------------------------------------------------------------------------------
# 1. 提取 Monocytes 亚群
# ------------------------------------------------------------------------------
cat("\n--- 步骤1: 提取 Monocytes 亚群 ---\n")
mono_cells <- subset(pbmc, subset = SingleR.labels == "Monocytes")
cat("Monocytes 细胞数:", ncol(mono_cells), "\n")
cat("Group 分布:\n")
print(table(mono_cells$Group))

# ------------------------------------------------------------------------------
# 2. 子集精细降维聚类
# ------------------------------------------------------------------------------
cat("\n--- 步骤2: 子集精细聚类 ---\n")
mono_cells <- JoinLayers(mono_cells)
mono_cells <- FindVariableFeatures(mono_cells, selection.method = "vst", nfeatures = 2000)
mono_cells <- ScaleData(mono_cells)
mono_cells <- RunPCA(mono_cells, verbose = FALSE)
mono_cells <- FindNeighbors(mono_cells, dims = 1:15)
mono_cells <- FindClusters(mono_cells,
                            resolution = subcluster_resolution,
                            algorithm = subcluster_algorithm)  # Leiden
mono_cells <- RunUMAP(mono_cells, dims = 1:15, verbose = FALSE)

cat("聚类方法:", subcluster_method, "| resolution:", subcluster_resolution, "\n")
cat("聚类数:", length(unique(Idents(mono_cells))), "\n")
cat("聚类分布:\n")
print(table(Idents(mono_cells)))

# 生成 mono_cluster_id：按实际 seurat_clusters 排序映射，兼容 Louvain(常从0开始) 和 Leiden(可从1开始)
cluster_ids <- sort(unique(as.character(mono_cells$seurat_clusters)))
cluster_map <- setNames(paste0("Mono_", seq_along(cluster_ids)), cluster_ids)
mono_cells$mono_cluster_id <- unname(cluster_map[as.character(mono_cells$seurat_clusters)])
mono_cells$mono_cluster_id <- factor(mono_cells$mono_cluster_id,
                                     levels = unname(cluster_map))
cat("cluster -> mono_cluster_id 映射:
")
print(data.frame(seurat_clusters = names(cluster_map), mono_cluster_id = unname(cluster_map)))

# ------------------------------------------------------------------------------
# 3. 计算 Marker 基因
# ------------------------------------------------------------------------------
cat("\n--- 步骤3: 计算 Monocyte 亚群 Marker ---\n")
Idents(mono_cells) <- "mono_cluster_id"

all_markers <- FindAllMarkers(mono_cells,
                              only.pos = TRUE,
                              min.pct = 0.25,
                              logfc.threshold = 0.25,
                              verbose = FALSE)

top30_markers <- all_markers %>%
  group_by(cluster) %>%
  slice_max(n = 30, order_by = avg_log2FC)

write.csv(all_markers,
          file.path(files_dir, paste0("introns_Markers_All_Monocytes_Subgroups_", subcluster_method, ".csv")),
          row.names = FALSE)
write.csv(top30_markers,
          file.path(files_dir, paste0("introns_Markers_Top30_Monocytes_Subgroups_", subcluster_method, ".csv")),
          row.names = FALSE)
cat("Marker 表已保存\n")

# ------------------------------------------------------------------------------
# 4. 局部坐标 UMAP (2x2 网格)
# ------------------------------------------------------------------------------
cat("\n--- 步骤4: 生成局部坐标 UMAP ---\n")

umap_coords <- Embeddings(mono_cells, reduction = "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])
x_buffer <- diff(x_range) * 0.05
y_buffer <- diff(y_range) * 0.05
global_xlim <- x_range + c(-x_buffer, x_buffer)
global_ylim <- y_range + c(-y_buffer, y_buffer)

p_total_local <- DimPlot(mono_cells, reduction = "umap", group.by = "mono_cluster_id",
                         label = TRUE, repel = TRUE) +
  ggtitle("Introns Monocyte Subgroups (Local Total)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_blank())

p_cold_local <- DimPlot(subset(mono_cells, subset = Group == "Cold_4C"), reduction = "umap",
                        group.by = "mono_cluster_id", label = FALSE) +
  ggtitle("Cold_4C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_rt_local <- DimPlot(subset(mono_cells, subset = Group == "RT_25C"), reduction = "umap",
                      group.by = "mono_cluster_id", label = FALSE) +
  ggtitle("RT_25C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_tn_local <- DimPlot(subset(mono_cells, subset = Group == "TN_30C"), reduction = "umap",
                      group.by = "mono_cluster_id", label = FALSE) +
  ggtitle("TN_30C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_final_local <- (p_total_local | p_cold_local) / (p_rt_local | p_tn_local) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.text = element_text(size = 9)) &
  coord_cartesian(xlim = global_xlim, ylim = global_ylim)

ggsave(file.path(plots_dir, paste0("introns_UMAP_Monocytes_Subgroups_Local_", subcluster_method, ".png")),
       plot = p_final_local, width = 15, height = 11, dpi = 300)
ggsave(file.path(plots_dir, paste0("introns_UMAP_Monocytes_Subgroups_Local_", subcluster_method, ".pdf")),
       plot = p_final_local, width = 15, height = 11)

# ------------------------------------------------------------------------------
# 5. 全局坐标 UMAP (映射回 PBMC)
# ------------------------------------------------------------------------------
cat("\n--- 步骤5: 生成全局坐标 UMAP ---\n")

pbmc$mono_viz_group <- "Others"
pbmc@meta.data[Cells(mono_cells), "mono_viz_group"] <- as.character(mono_cells$mono_cluster_id)
mono_levels <- sort(unique(as.character(mono_cells$mono_cluster_id)))
pbmc$mono_viz_group <- factor(pbmc$mono_viz_group, levels = c(mono_levels, "Others"))

cells_cold <- Cells(subset(mono_cells, subset = Group == "Cold_4C"))
cells_rt <- Cells(subset(mono_cells, subset = Group == "RT_25C"))
cells_tn <- Cells(subset(mono_cells, subset = Group == "TN_30C"))

p_total_cls <- DimPlot(pbmc, cells = Cells(mono_cells), reduction = "umap",
                       group.by = "mono_viz_group", label = TRUE, repel = TRUE) +
  ggtitle("Introns Monocyte Subgroups (Global Total)") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_blank())

p_cold_cls <- DimPlot(pbmc, cells = cells_cold, reduction = "umap",
                      group.by = "mono_viz_group", label = FALSE) +
  ggtitle("Cold_4C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_rt_cls <- DimPlot(pbmc, cells = cells_rt, reduction = "umap",
                    group.by = "mono_viz_group", label = FALSE) +
  ggtitle("RT_25C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_tn_cls <- DimPlot(pbmc, cells = cells_tn, reduction = "umap",
                    group.by = "mono_viz_group", label = FALSE) +
  ggtitle("TN_30C") +
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
  NoLegend()

p_final_cls <- (p_total_cls | p_cold_cls) / (p_rt_cls | p_tn_cls) +
  plot_layout(guides = "collect", axes = "collect") &
  theme(legend.text = element_text(size = 9))

ggsave(file.path(plots_dir, paste0("introns_UMAP_Monocytes_Subgroups_Global_", subcluster_method, ".png")),
       plot = p_final_cls, width = 15, height = 11, dpi = 300)

# 保存 monocyte 子集对象
qsave(mono_cells, file.path(objects_dir, paste0("coldmouse_introns_pbmc_monocytes_subclustered_", subcluster_method, ".qs")))

# ------------------------------------------------------------------------------
# 6. Monocytes 组间 DEG: Cold vs RT
# ------------------------------------------------------------------------------
cat("\n--- 步骤6: Cold vs RT DEG ---\n")
Idents(mono_cells) <- "Group"

deg_cold_vs_rt <- FindMarkers(mono_cells,
                              ident.1 = "Cold_4C", ident.2 = "RT_25C",
                              logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
deg_cold_vs_rt <- deg_cold_vs_rt %>% rownames_to_column(var = "gene") %>% arrange(desc(avg_log2FC))
write.csv(deg_cold_vs_rt,
          file.path(files_dir, paste0("introns_DEG_Monocytes_Cold_vs_RT_", subcluster_method, ".csv")),
          row.names = FALSE)
cat("Cold vs RT DEG 完成 (", nrow(deg_cold_vs_rt), "个基因 )\n")

# ------------------------------------------------------------------------------
# 7. Monocytes 组间 DEG: Cold vs TN
# ------------------------------------------------------------------------------
cat("\n--- 步骤7: Cold vs TN DEG ---\n")
deg_cold_vs_tn <- FindMarkers(mono_cells,
                              ident.1 = "Cold_4C", ident.2 = "TN_30C",
                              logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
deg_cold_vs_tn <- deg_cold_vs_tn %>% rownames_to_column(var = "gene") %>% arrange(desc(avg_log2FC))
write.csv(deg_cold_vs_tn,
          file.path(files_dir, paste0("introns_DEG_Monocytes_Cold_vs_TN_", subcluster_method, ".csv")),
          row.names = FALSE)
cat("Cold vs TN DEG 完成 (", nrow(deg_cold_vs_tn), "个基因 )\n")

# ------------------------------------------------------------------------------
# 8. 各 cluster vs Others DEG
# ------------------------------------------------------------------------------
cat("\n--- 步骤8: 各 cluster vs Others DEG ---\n")
Idents(mono_cells) <- "mono_cluster_id"

for (cl in levels(mono_cells$mono_cluster_id)) {
  cat("  计算", cl, "vs Others...\n")
  deg_tmp <- FindMarkers(mono_cells, ident.1 = cl, ident.2 = NULL,
                         logfc.threshold = 0.25, min.pct = 0.1, verbose = FALSE)
  deg_tmp <- deg_tmp %>%
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC))

  write.csv(deg_tmp,
            file.path(files_dir, paste0("introns_DEG_", cl, "_vs_Others_", subcluster_method, ".csv")),
            row.names = FALSE)

  top20 <- deg_tmp %>% slice_max(n = 20, order_by = avg_log2FC, with_ties = FALSE)
  write.csv(top20,
            file.path(files_dir, paste0("introns_DEG_", cl, "_vs_Others_Top20_", subcluster_method, ".csv")),
            row.names = FALSE)
}

cat("\n=== PBMC Monocyte 亚群分析完成 ===\n")
cat("对象保存:", file.path(objects_dir, paste0("coldmouse_introns_pbmc_monocytes_subclustered_", subcluster_method, ".qs")), "\n")
cat("结果目录:", files_dir, "\n")
cat("图片目录:", plots_dir, "\n")
