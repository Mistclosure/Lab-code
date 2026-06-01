# ==============================================================================
# utils_qc.R: 质量控制工具
# ==============================================================================
# 提供 QC 计算、过滤、核糖体基因处理、双细胞检测等功能。
# ==============================================================================

#' 计算 QC 指标
#'
#' @param obj Seurat 对象
#' @param ribo_pattern 核糖体基因模式
#' @return 添加了 QC 指标的 Seurat 对象
calculate_qc_metrics <- function(obj, ribo_pattern = "^Rp[sl]") {
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = ribo_pattern)
  obj
}

#' 应用 QC 过滤
#'
#' @param obj Seurat 对象
#' @param min_features 最少基因数
#' @param max_features 最多基因数 (NULL 表示不设上限)
#' @param max_count 最大 UMI 数 (NULL 表示不设上限)
#' @param percent_mt 线粒体基因百分比阈值
#' @return 过滤后的 Seurat 对象
apply_qc_filter <- function(obj, min_features = 200, max_features = NULL,
                            max_count = NULL, percent_mt = 15) {
  subset_str <- paste0("nFeature_RNA >= ", min_features)
  if (!is.null(max_features)) {
    subset_str <- paste0(subset_str, " & nFeature_RNA <= ", max_features)
  }
  if (!is.null(max_count)) {
    subset_str <- paste0(subset_str, " & nCount_RNA <= ", max_count)
  }
  subset_str <- paste0(subset_str, " & percent.mt < ", percent_mt)

  subset(obj, subset = eval(parse(text = subset_str)))
}

#' 计算 QC 分位数表
#'
#' @param obj Seurat 对象
#' @return QC 分位数 data.frame
calculate_qc_quantiles <- function(obj) {
  data.frame(
    metric = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
    min = c(
      min(obj$nFeature_RNA),
      min(obj$nCount_RNA),
      min(obj$percent.mt),
      min(obj$percent.rb)
    ),
    q25 = c(
      quantile(obj$nFeature_RNA, 0.25),
      quantile(obj$nCount_RNA, 0.25),
      quantile(obj$percent.mt, 0.25),
      quantile(obj$percent.rb, 0.25)
    ),
    median = c(
      median(obj$nFeature_RNA),
      median(obj$nCount_RNA),
      median(obj$percent.mt),
      median(obj$percent.rb)
    ),
    q75 = c(
      quantile(obj$nFeature_RNA, 0.75),
      quantile(obj$nCount_RNA, 0.75),
      quantile(obj$percent.mt, 0.75),
      quantile(obj$percent.rb, 0.75)
    ),
    max = c(
      max(obj$nFeature_RNA),
      max(obj$nCount_RNA),
      max(obj$percent.mt),
      max(obj$percent.rb)
    ),
    stringsAsFactors = FALSE
  )
}

#' 移除核糖体基因
#'
#' @param obj Seurat 对象
#' @param ribo_pattern 核糖体基因模式
#' @return 移除了核糖体基因的 Seurat 对象
remove_ribo_genes <- function(obj, ribo_pattern = "^Rp[sl]") {
  ribo_genes <- grep(ribo_pattern, rownames(obj), value = TRUE, ignore.case = TRUE)
  non_ribo_genes <- rownames(obj)[!rownames(obj) %in% ribo_genes]
  subset(obj, features = non_ribo_genes)
}

#' 运行 scDblFinder 去双胞
#'
#' @param obj Seurat 对象
#' @param resolution 聚类分辨率
#' @param pca_dims PCA 维度
#' @return 去双胞后的 Seurat 对象
run_scDblFinder <- function(obj, resolution = 0.5, pca_dims = 20) {
  if (ncol(obj) < 50) {
    message("细胞数 < 50，跳过 scDblFinder")
    return(obj)
  }

  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:pca_dims)
  obj <- FindClusters(obj, resolution = resolution)

  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class

    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    message("   [scDblFinder] 双细胞数: ", n_dbl, " (", round(n_dbl / ncol(obj) * 100, 2), "%)")

    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    message("   [Filter] 去双胞后细胞数: ", ncol(obj))
  }, error = function(e) {
    message("   scDblFinder 警告: ", e$message, " - 保留所有细胞继续")
  })

  obj
}

#' 生成 QC 小提琴图和散点图
#'
#' @param obj Seurat 对象
#' @param output_dir 输出目录
#' @param prefix 文件名前缀
#' @param width 图片宽度
#' @param height 图片高度
#' @param dpi 图片分辨率
generate_qc_plots <- function(obj, output_dir, prefix = "coldmouse_introns",
                              width = 12, height = 8, dpi = 300) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 小提琴图
  p_vln <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
                    ncol = 4, pt.size = 0)
  ggsave(file.path(output_dir, paste0(prefix, "_qc_violin.png")),
         plot = p_vln, width = width, height = height, dpi = dpi)

  # 散点图
  p_scatter1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  p_scatter2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  p_scatter <- p_scatter1 + p_scatter2
  ggsave(file.path(output_dir, paste0(prefix, "_qc_scatter.png")),
         plot = p_scatter, width = width, height = height / 2, dpi = dpi)

  message("QC 图片已保存到: ", output_dir)
}
