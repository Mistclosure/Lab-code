# ==============================================================================
# utils_seurat.R: Seurat 流程工具
# ==============================================================================
# 提供标准化、降维、聚类、注释等 Seurat 流程功能。
# ==============================================================================

#' 运行标准 Seurat 流程 (无 Harmony)
#'
#' @param obj Seurat 对象
#' @param tissue_name 组织名称
#' @param method_id 聚类算法 ID (1=Louvain, 4=Leiden)
#' @param n_pcs PCA 维度
#' @param resolution 聚类分辨率
#' @param nfeatures 高变基因数量
#' @return 处理后的 Seurat 对象
run_seurat_standard <- function(obj, tissue_name, method_id = 1,
                                n_pcs = 50, resolution = 1.0, nfeatures = 2000) {
  message(">>> 处理组织: ", tissue_name, " | 细胞数: ", ncol(obj))
  if (ncol(obj) < 50) {
    message("   细胞数 < 50，跳过")
    return(NULL)
  }

  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:n_pcs)
  obj <- FindClusters(obj, resolution = resolution, algorithm = method_id)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:n_pcs)

  obj
}

#' 运行带 Harmony 的 Seurat 流程
#'
#' @param obj Seurat 对象
#' @param tissue_name 组织名称
#' @param method_id 聚类算法 ID
#' @param n_pcs PCA 维度
#' @param resolution 聚类分辨率
#' @param harmony_group Harmony 分组变量
#' @return 处理后的 Seurat 对象
run_seurat_harmony <- function(obj, tissue_name, method_id = 1,
                               n_pcs = 50, resolution = 1.0, harmony_group = "orig.ident") {
  message(">>> 处理组织: ", tissue_name, " (Harmony) | 细胞数: ", ncol(obj))
  if (ncol(obj) < 50) {
    message("   细胞数 < 50，跳过")
    return(NULL)
  }

  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = n_pcs, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = harmony_group, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:n_pcs)
  obj <- FindClusters(obj, resolution = resolution, algorithm = method_id)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:n_pcs)

  obj
}

#' 运行 SingleR 全局注释
#'
#' @param obj Seurat 对象
#' @param ref_data SingleR 参考数据集
#' @param labels_col 参考数据的标签列
#' @return 添加了 SingleR.labels 的 Seurat 对象
run_singler_annotation <- function(obj, ref_data, labels_col = "label.main") {
  obj <- JoinLayers(obj)
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  pred_res <- SingleR(test = expr_mat, ref = ref_data, labels = ref_data[[labels_col]])
  obj$SingleR.labels <- pred_res$labels
  obj
}

#' 计算并导出 marker 表
#'
#' @param obj Seurat 对象
#' @param tissue_name 组织名称
#' @param cluster_method 聚类方法名
#' @param output_dir 输出目录
#' @param min_cells_threshold 最小细胞数阈值
#' @param prefix 文件名前缀
#' @return marker 表 data.frame
export_markers <- function(obj, tissue_name, cluster_method = "louvain",
                           output_dir, min_cells_threshold = 20,
                           prefix = "Introns") {
  # 过滤小亚群
  cell_counts <- table(obj$SingleR.labels)
  keep_labels <- names(cell_counts[cell_counts >= min_cells_threshold])
  if (length(keep_labels) == 0) {
    message("   ", tissue_name, " 中没有符合条件的亚群，跳过")
    return(NULL)
  }
  obj <- subset(obj, subset = SingleR.labels %in% keep_labels)

  # 计算 marker
  Idents(obj) <- "seurat_clusters"
  all_markers <- FindAllMarkers(obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                verbose = FALSE)

  # 建立 Cluster 与 SingleR 标签的映射
  cluster_to_label <- as.data.frame(obj@meta.data) %>%
    select(seurat_clusters, SingleR.labels) %>%
    mutate(SingleR.labels = as.character(unlist(SingleR.labels))) %>%
    group_by(seurat_clusters) %>%
    dplyr::count(SingleR.labels) %>%
    slice_max(n = 1, order_by = n, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster = seurat_clusters, SingleR_Annotation = SingleR.labels)

  all_markers <- all_markers %>%
    left_join(cluster_to_label, by = "cluster") %>%
    mutate(cluster_with_anno = paste0("Cluster ", cluster, " (", SingleR_Annotation, ")"))

  top20_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)

  # 保存
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(all_markers,
            file.path(output_dir, paste0(prefix, "_Markers_All_", cluster_method, "_", tissue_name, ".csv")),
            row.names = FALSE)
  write.csv(top20_markers,
            file.path(output_dir, paste0(prefix, "_Markers_Top20_", cluster_method, "_", tissue_name, ".csv")),
            row.names = FALSE)

  message("   Marker 表已导出: ", tissue_name)
  all_markers
}
