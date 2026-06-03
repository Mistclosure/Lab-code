# ==============================================================================
# 01_export_malignant_qs_to_mtx.R
# 从 Seurat qs 对象导出 counts matrix + metadata + reductions
# 不使用 reticulate / zellkonverter
# ==============================================================================

suppressPackageStartupMessages({
  library(qs)
  library(Seurat)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# 1. 路径设置
# ------------------------------------------------------------------------------

work_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"

input_qs <- file.path(work_dir, "Malignant_RNA_assay.qs")

out_dir <- file.path(work_dir, "Malignant_RNA_assay_export_for_python")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("输入 qs 文件：\n")
cat(input_qs, "\n\n")

cat("输出目录：\n")
cat(out_dir, "\n\n")

# ------------------------------------------------------------------------------
# 2. 读取 Seurat 对象
# ------------------------------------------------------------------------------

cat("正在读取 Seurat qs 对象...\n")

obj <- qread(input_qs)

cat("读取完成。\n")
cat("对象类型：\n")
print(class(obj))

if (!inherits(obj, "Seurat")) {
  stop("读取出来的对象不是 Seurat object，请检查 Malignant_RNA_assay.qs。")
}

cat("原始细胞数：", ncol(obj), "\n")
cat("原始基因数：", nrow(obj), "\n\n")

# ------------------------------------------------------------------------------
# 3. 设置 RNA assay
# ------------------------------------------------------------------------------

assay_names <- names(obj@assays)

cat("当前 assay：\n")
print(assay_names)

if (!"RNA" %in% assay_names) {
  stop("对象中没有 RNA assay。当前 assay 为：", paste(assay_names, collapse = ", "))
}

DefaultAssay(obj) <- "RNA"

# Seurat v5 如果存在多个 layer，先 JoinLayers
cat("\n尝试 JoinLayers...\n")

obj <- tryCatch(
  {
    JoinLayers(obj)
  },
  error = function(e) {
    cat("JoinLayers 跳过：", conditionMessage(e), "\n")
    obj
  }
)

# ------------------------------------------------------------------------------
# 4. 提取 counts matrix
# ------------------------------------------------------------------------------

get_assay_matrix <- function(seurat_obj, assay = "RNA", layer_or_slot = "counts") {
  mat <- tryCatch(
    {
      GetAssayData(seurat_obj, assay = assay, layer = layer_or_slot)
    },
    error = function(e1) {
      GetAssayData(seurat_obj, assay = assay, slot = layer_or_slot)
    }
  )
  return(mat)
}

cat("\n正在提取 RNA raw counts...\n")

counts <- get_assay_matrix(
  seurat_obj = obj,
  assay = "RNA",
  layer_or_slot = "counts"
)

if (is.null(counts) || nrow(counts) == 0 || ncol(counts) == 0) {
  stop("counts matrix 为空，无法继续。")
}

counts <- as(counts, "dgCMatrix")

cat("counts 原始维度 genes × cells：\n")
print(dim(counts))

# ------------------------------------------------------------------------------
# 5. 过滤全零基因和全零细胞
# ------------------------------------------------------------------------------

cat("\n正在过滤全零基因和全零细胞...\n")

keep_genes <- Matrix::rowSums(counts) > 0
keep_cells <- Matrix::colSums(counts) > 0

counts <- counts[keep_genes, keep_cells, drop = FALSE]

cat("counts 过滤后维度 genes × cells：\n")
print(dim(counts))

# ------------------------------------------------------------------------------
# 6. 导出 counts matrix
# ------------------------------------------------------------------------------

cat("\n正在导出 counts matrix 为 Matrix Market 格式...\n")

Matrix::writeMM(
  obj = counts,
  file = file.path(out_dir, "counts_genes_by_cells.mtx")
)

# ------------------------------------------------------------------------------
# 7. 导出 genes.tsv 和 barcodes.tsv
# ------------------------------------------------------------------------------

cat("正在导出 genes.tsv 和 barcodes.tsv...\n")

gene_symbol <- rownames(counts)
gene_id <- make.unique(gene_symbol)

genes_df <- data.frame(
  gene_id = gene_id,
  gene_symbol = gene_symbol,
  stringsAsFactors = FALSE
)

barcodes_df <- data.frame(
  cell = colnames(counts),
  stringsAsFactors = FALSE
)

write.table(
  genes_df,
  file = file.path(out_dir, "genes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

write.table(
  barcodes_df,
  file = file.path(out_dir, "barcodes.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ------------------------------------------------------------------------------
# 8. 导出 metadata
# ------------------------------------------------------------------------------

cat("正在导出 metadata.tsv...\n")

meta <- obj@meta.data[colnames(counts), , drop = FALSE]
meta <- as.data.frame(meta, check.names = FALSE)

clean_meta_column <- function(x) {
  if (is.factor(x)) {
    return(as.character(x))
  }
  if (inherits(x, c("Date", "POSIXct", "POSIXlt"))) {
    return(as.character(x))
  }
  if (is.list(x)) {
    return(vapply(x, function(y) paste(as.character(y), collapse = ";"), character(1)))
  }
  return(x)
}

meta[] <- lapply(meta, clean_meta_column)

meta$cell <- rownames(meta)
meta <- meta[, c("cell", setdiff(colnames(meta), "cell")), drop = FALSE]

write.table(
  meta,
  file = file.path(out_dir, "metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# ------------------------------------------------------------------------------
# 9. 可选：导出 Seurat 的 normalized data
# 默认关闭，因为 cNMF 推荐使用 raw counts
# ------------------------------------------------------------------------------

EXPORT_LOGDATA <- FALSE

if (EXPORT_LOGDATA) {
  cat("\n正在尝试导出 RNA data matrix...\n")
  
  data_mat <- tryCatch(
    {
      get_assay_matrix(obj, assay = "RNA", layer_or_slot = "data")
    },
    error = function(e) {
      NULL
    }
  )
  
  if (!is.null(data_mat)) {
    data_mat <- as(data_mat, "dgCMatrix")
    data_mat <- data_mat[gene_symbol, colnames(counts), drop = FALSE]
    
    Matrix::writeMM(
      obj = data_mat,
      file = file.path(out_dir, "logdata_genes_by_cells.mtx")
    )
    
    cat("logdata 已导出。\n")
  } else {
    cat("未找到 RNA data matrix，跳过。\n")
  }
}

# ------------------------------------------------------------------------------
# 10. 导出降维结果：pca / harmony / umap / tsne 等
# ------------------------------------------------------------------------------

cat("\n正在导出 reductions...\n")

reduction_names <- names(obj@reductions)

cat("检测到 reductions：\n")
print(reduction_names)

for (red in reduction_names) {
  
  emb <- tryCatch(
    {
      Embeddings(obj, reduction = red)
    },
    error = function(e) {
      NULL
    }
  )
  
  if (is.null(emb)) {
    next
  }
  
  if (!all(colnames(counts) %in% rownames(emb))) {
    cat("跳过 reduction：", red, "，因为细胞名不完全匹配。\n")
    next
  }
  
  emb <- emb[colnames(counts), , drop = FALSE]
  
  emb_df <- as.data.frame(emb, check.names = FALSE)
  emb_df$cell <- rownames(emb_df)
  emb_df <- emb_df[, c("cell", setdiff(colnames(emb_df), "cell")), drop = FALSE]
  
  out_file <- file.path(out_dir, paste0("reduction_", red, ".tsv"))
  
  write.table(
    emb_df,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )
  
  cat("已导出 reduction：", red, "\n")
}

# ------------------------------------------------------------------------------
# 11. 输出 summary
# ------------------------------------------------------------------------------

cat("\n================ R 导出完成 ================\n")
cat("输出目录：", out_dir, "\n")
cat("counts 文件：counts_genes_by_cells.mtx\n")
cat("genes 文件：genes.tsv\n")
cat("barcodes 文件：barcodes.tsv\n")
cat("metadata 文件：metadata.tsv\n")
cat("最终细胞数：", ncol(counts), "\n")
cat("最终基因数：", nrow(counts), "\n")
cat("============================================\n")