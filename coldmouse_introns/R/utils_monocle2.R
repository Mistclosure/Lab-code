# ==============================================================================
# utils_monocle2.R: Monocle 2 轨迹分析工具
# ==============================================================================
# 提供 Seurat 到 Monocle 2 的转换、轨迹推断等功能。
# 注意：monocle2 包不是 monocle3。
# ==============================================================================

#' 将 Seurat 对象转换为 Monocle 2 CellDataSet
#'
#' @param obj Seurat 对象
#' @param assay 使用的 assay 名称
#' @return Monocle 2 CellDataSet 对象
seurat_to_monocle2 <- function(obj, assay = "RNA") {
  counts <- GetAssayData(obj, assay = assay, layer = "counts")
  pd <- new("AnnotatedDataFrame", data = obj@meta.data)
  fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  fd <- new("AnnotatedDataFrame", data = fd)

  cds <- newCellDataSet(
    as(counts, "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    expressionFamily = negbinomial.size()
  )
  cds
}

#' 运行 Monocle 2 标准轨迹分析流程
#'
#' @param cds Monocle 2 CellDataSet
#' @param root_marker 根节点 marker 基因
#' @param num_dim 降维维度
#' @return 处理后的 CellDataSet
run_monocle2_trajectory <- function(cds, root_marker = "Ly6c2", num_dim = 50) {
  # 过滤低表达基因
  cds <- estimateSizeFactors(cds)
  cds <- estimateDispersions(cds)

  # 检测差异表达基因
  cds <- detectGenes(cds, min_expr = 0.1)
  expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

  # 找差异基因
  diff_test_res <- differentialGeneTest(cds[expressed_genes, ],
                                         fullModelFormulaModel = "~Tissue")
  ordering_genes <- rownames(diff_test_res[order(diff_test_res$qval), ])[1:min(1000, nrow(diff_test_res))]

  # 设置排序基因
  cds <- setOrderingFilter(cds, ordering_genes)

  # 降维
  cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")

  # 排序
  if (root_marker %in% rownames(cds)) {
    root_cells <- which.max(exprs(cds)[root_marker, ])
    cds <- orderCells(cds, root_state = root_cells)
  } else {
    cds <- orderCells(cds)
  }

  cds
}

#' 生成 Monocle 2 轨迹图
#'
#' @param cds Monocle 2 CellDataSet
#' @param color_by 颜色分组变量
#' @param output_dir 输出目录
#' @param prefix 文件名前缀
#' @param width 图片宽度
#' @param height 图片高度
#' @param dpi 图片分辨率
generate_monocle2_plots <- function(cds, color_by = "State",
                                    output_dir, prefix = "introns",
                                    width = 8, height = 6, dpi = 300) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # 轨迹图
  p <- plot_cell_trajectory(cds, color_by = color_by) +
    ggtitle(paste(prefix, "Trajectory"))
  ggsave(file.path(output_dir, paste0(prefix, "_trajectory_", color_by, ".png")),
         plot = p, width = width, height = height, dpi = dpi)

  # 伪时间图
  p_time <- plot_cell_trajectory(cds, color_by = "Pseudotime") +
    ggtitle(paste(prefix, "Pseudotime"))
  ggsave(file.path(output_dir, paste0(prefix, "_pseudotime.png")),
         plot = p_time, width = width, height = height, dpi = dpi)

  message("Monocle 2 轨迹图已保存到: ", output_dir)
}
