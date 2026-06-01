# ==============================================================================
# utils_io.R: 输入输出工具
# ==============================================================================
# 提供 loom 文件读取、Seurat 对象保存/加载、日志记录等功能。
# ==============================================================================

#' 读取 introns loom 文件并返回稀疏矩阵
#'
#' @param loom_file loom 文件路径
#' @param sample 样本 ID
#' @param chunk_cells 每次读取的细胞数
#' @return 稀疏矩阵 (gene x cell)
read_introns_loom <- function(loom_file, sample, chunk_cells = 5000) {
  if (!file.exists(loom_file)) {
    stop("Loom 文件未找到: ", loom_file)
  }

  f <- H5File$new(loom_file, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  ds <- f[["matrix"]]
  dims <- ds$dims
  n_cells <- dims[1]
  n_genes <- dims[2]

  genes <- as.character(f[["row_attrs/Gene"]][])
  barcodes <- as.character(f[["col_attrs/obs_names"]][])

  if (length(genes) != n_genes) {
    stop("基因数量与矩阵列数不匹配: ", loom_file)
  }
  if (length(barcodes) != n_cells) {
    stop("细胞 barcode 数量与矩阵行数不匹配: ", loom_file)
  }

  genes <- make.unique(genes)
  barcodes <- gsub("_", "-", barcodes)
  barcodes <- sub(paste0("^", sample, ":"), "", barcodes)

  counts_chunks <- vector("list", ceiling(n_cells / chunk_cells))
  chunk_id <- 1
  total_chunks <- length(counts_chunks)

  for (start_cell in seq(1, n_cells, by = chunk_cells)) {
    end_cell <- min(start_cell + chunk_cells - 1, n_cells)
    cell_idx <- start_cell:end_cell
    message(
      "   读取 ", sample, " loom 分块 ", chunk_id, "/", total_chunks,
      "：细胞 ", start_cell, "-", end_cell
    )

    block <- ds[cell_idx, ]
    nz <- which(block > 0, arr.ind = TRUE)

    if (nrow(nz) > 0) {
      counts_chunks[[chunk_id]] <- sparseMatrix(
        i = nz[, 2],
        j = nz[, 1],
        x = as.numeric(block[nz]),
        dims = c(n_genes, length(cell_idx))
      )
    } else {
      counts_chunks[[chunk_id]] <- sparseMatrix(
        i = integer(0),
        j = integer(0),
        x = numeric(0),
        dims = c(n_genes, length(cell_idx))
      )
    }

    rm(block, nz)
    gc(verbose = FALSE)
    chunk_id <- chunk_id + 1
  }

  counts <- do.call(cbind, counts_chunks)
  rownames(counts) <- genes
  colnames(counts) <- barcodes
  counts
}

#' 检查 loom 文件结构
#'
#' @param loom_file loom 文件路径
#' @return 包含结构信息的列表
check_loom_structure <- function(loom_file) {
  if (!file.exists(loom_file)) {
    stop("Loom 文件未找到: ", loom_file)
  }
  f <- H5File$new(loom_file, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  ds <- f[["matrix"]]
  dims <- ds$dims

  genes <- as.character(f[["row_attrs/Gene"]][])
  barcodes <- as.character(f[["col_attrs/obs_names"]][])

  list(
    n_cells = dims[1],
    n_genes = dims[2],
    genes = genes,
    barcodes = barcodes,
    matrix_dims = dims,
    is_gene_x_cell = dims[2] == length(genes)
  )
}

#' 保存 Seurat 对象
#'
#' @param obj Seurat 对象
#' @param file_path 保存路径
#' @param overwrite 是否覆盖已存在的文件
save_seurat <- function(obj, file_path, overwrite = FALSE) {
  if (file.exists(file_path) && !overwrite) {
    warning("文件已存在，跳过保存: ", file_path)
    return(invisible(NULL))
  }
  qsave(obj, file_path)
  message("已保存: ", file_path)
}

#' 加载 Seurat 对象
#'
#' @param file_path 文件路径
#' @return Seurat 对象
load_seurat <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("文件未找到: ", file_path)
  }
  qread(file_path)
}
