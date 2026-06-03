# ==============================================================================
# 03_clustering_and_global_anno_harmony_introns.R
# coldmouse_introns：带 Harmony 的降维聚类、SingleR 全局注释和 marker 输出
# ==============================================================================
# 输入：coldmouse_introns_sc_combined_qc_Cleaned.qs
# 输出：coldmouse_introns_sc_by_tissue_harmony_louvain.qs、marker csv、UMAP png
# 参考：docs/original_reference/clustering_and_global_anno.r
# 注意：默认不运行，只有 config run_harmony=true 时使用
# ==============================================================================

# 加载配置和工具
PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))
source(file.path(PROJECT_HOME, "R", "utils_seurat.R"))

config <- load_config()
setup_performance(config)

# 检查是否启用 Harmony
if (!config$harmony$run_harmony) {
  message("Harmony 未启用 (config harmony.run_harmony = false)。跳过此脚本。")
  message("如需运行 Harmony，请在 config.yaml 中设置 harmony.run_harmony = true")
  quit(save = "no", status = 0)
}

# 加载 R 包
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(leidenbase)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(SingleR)
library(celldex)

# 设置输出目录
objects_dir <- get_results_dir(config, "objects")
files_dir <- get_results_dir(config, "files")
plots_dir <- get_results_dir(config, "plots")

# 聚类参数
cluster_config <- config$clustering
cluster_method <- cluster_config$method
algo_id <- ifelse(cluster_method == "leiden", 4, 1)
resolution <- cluster_config$resolution
n_pcs <- cluster_config$pca_dims
nfeatures <- cluster_config$nfeatures

# Harmony 参数
harmony_config <- config$harmony
harmony_group <- harmony_config$group_by_vars

# SingleR 参数
singler_config <- config$singler
min_cells_threshold <- singler_config$min_cells_threshold

# 可视化参数
viz_config <- config$visualization

# 读取 QC 后数据
qc_file <- file.path(objects_dir, "coldmouse_introns_sc_combined_qc_Cleaned.qs")
if (!file.exists(qc_file)) {
  stop("未找到 introns 质控输出文件：", qc_file)
}
sc_combined <- load_seurat(qc_file)

# ------------------------------------------------------------------------------
# 1. 按组织独立进行标准化、降维和聚类，使用 Harmony
# ------------------------------------------------------------------------------
print(paste0("步骤1/4：使用 ", cluster_method, " + Harmony 进行降维与聚类。"))
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

sc_by_tissue <- lapply(names(sc_by_tissue), function(x) {
  run_seurat_harmony(sc_by_tissue[[x]], x,
                     method_id = algo_id,
                     n_pcs = n_pcs,
                     resolution = resolution,
                     harmony_group = harmony_group)
})
names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow")
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]

if (length(sc_by_tissue) == 0) {
  stop("没有可用于后续分析的组织对象。")
}

# ------------------------------------------------------------------------------
# 2. 运行 SingleR 对所有细胞做全局注释
# ------------------------------------------------------------------------------
print("步骤2/4：运行 SingleR 对所有细胞进行全局注释。")
mouse_ref <- celldex::MouseRNAseqData()

for (i in seq_along(sc_by_tissue)) {
  tissue_name <- names(sc_by_tissue)[i]
  print(paste(">>> 正在运行 SingleR 注释:", tissue_name))
  sc_by_tissue[[i]] <- run_singler_annotation(sc_by_tissue[[i]], mouse_ref)
  print(paste("   完成:", tissue_name, "全局注释"))
}

# 保存按组织对象
sc_by_tissue_file <- file.path(objects_dir,
                               paste0("coldmouse_introns_sc_by_tissue_harmony_", cluster_method, ".qs"))
save_seurat(sc_by_tissue, sc_by_tissue_file)

# ------------------------------------------------------------------------------
# 3. 绘制 UMAP 并导出 marker 表
# ------------------------------------------------------------------------------
print("步骤3/4：绘制 UMAP 并导出 marker 表。")

for (tissue_name in names(sc_by_tissue)) {
  print(paste(">>> 正在处理组织:", tissue_name))
  obj <- sc_by_tissue[[tissue_name]]

  # 导出 marker 表
  export_markers(obj, tissue_name,
                 cluster_method = paste0(cluster_method, "_harmony"),
                 output_dir = files_dir,
                 min_cells_threshold = min_cells_threshold,
                 prefix = "Introns")

  # 固定 SingleR 标签顺序
  all_labels <- sort(unique(obj$SingleR.labels))
  obj$SingleR.labels <- factor(obj$SingleR.labels, levels = all_labels)
  sc_by_tissue[[tissue_name]] <- obj

  # 绘制 UMAP (按 Annotation)
  p_total_anno <- DimPlot(obj, reduction = "umap", group.by = "SingleR.labels",
                          label = TRUE, repel = TRUE) +
    ggtitle(paste(tissue_name, "-", cluster_method, "(Introns Harmony Annotation)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_blank())

  p_cold_anno <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap",
                         group.by = "SingleR.labels", label = FALSE) +
    ggtitle("Cold_4C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_rt_anno <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap",
                       group.by = "SingleR.labels", label = FALSE) +
    ggtitle("RT_25C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_tn_anno <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap",
                       group.by = "SingleR.labels", label = FALSE) +
    ggtitle("TN_30C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.text = element_text(size = 9))

  ggsave(file.path(plots_dir, paste0("Introns_UMAP_Grid_Annotation_", cluster_method, "_harmony_", tissue_name, ".png")),
         plot = p_final_anno,
         width = viz_config$umap_width, height = viz_config$umap_height, dpi = viz_config$umap_dpi)

  # 绘制 UMAP (按 Clusters)
  p_total_cls <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters",
                         label = TRUE, repel = TRUE) +
    ggtitle(paste(tissue_name, "-", cluster_method, "(Introns Harmony Clusters)")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_blank())

  p_cold_cls <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap",
                        group.by = "seurat_clusters", label = FALSE) +
    ggtitle("Cold_4C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_rt_cls <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap",
                      group.by = "seurat_clusters", label = FALSE) +
    ggtitle("RT_25C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_tn_cls <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap",
                      group.by = "seurat_clusters", label = FALSE) +
    ggtitle("TN_30C") +
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()

  p_final_cls <- (p_total_cls | p_cold_cls) / (p_rt_cls | p_tn_cls) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.text = element_text(size = 9))

  ggsave(file.path(plots_dir, paste0("Introns_UMAP_Grid_Clusters_", cluster_method, "_harmony_", tissue_name, ".png")),
         plot = p_final_cls,
         width = viz_config$umap_width, height = viz_config$umap_height, dpi = viz_config$umap_dpi)

  print(paste("   完成:", tissue_name, "Annotation 和 Clusters 两版 UMAP 已保存。"))
}

# ------------------------------------------------------------------------------
# 4. 保存带筛选后对象的 introns 版本组织列表
# ------------------------------------------------------------------------------
save_seurat(sc_by_tissue,
            file.path(objects_dir,
                      paste0("coldmouse_introns_sc_by_tissue_harmony_", cluster_method, "_with_markers_input.qs")))
print("步骤4/4：脚本结束。请在 results/ 中查看 introns 版本 Harmony 输出。")
