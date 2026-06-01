# ==============================================================================
# coldmouse_introns：数据读取、质量控制、核糖体基因移除和双细胞过滤
# ==============================================================================
#
# 输入路径：
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns/rawdata
#
# 输出路径：
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns
#
# 参考脚本：
#   /home/zerui/code/coldmouse/coldmouse_analysis/using/general_preprocessing/data_loading_and_qc.r
#
# 主要修改点：
#   1. 将原始 10X matrix 文件夹读取方式改为读取带 introns 信息的 loom 文件。
#   2. 保留原项目的样本命名、metadata 映射、QC 阈值、核糖体基因移除、
#      Seurat v5 JoinLayers 和 scDblFinder 去双胞逻辑。
#   3. 输出文件名加入 coldmouse_introns/introns 标识，避免覆盖原 coldmouse 项目结果。
#   4. 本脚本只负责质控和预处理，不运行后续聚类、注释或 marker 分析。
# ==============================================================================

setwd("/mnt/disk1/qiuzerui/expriments/coldmouse_introns")

dir.create("files", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 0. 加载 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(patchwork)
library(hdf5r)

if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

# ------------------------------------------------------------------------------
# 1. 路径和 metadata 映射
# ------------------------------------------------------------------------------
project_dir <- "/mnt/disk1/qiuzerui/expriments/coldmouse_introns"
data_dir <- file.path(project_dir, "rawdata")
output_file <- file.path(project_dir, "coldmouse_introns_sc_combined_qc_Cleaned.qs")

sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")
tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")
group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")

# ------------------------------------------------------------------------------
# 2. loom 文件读取和 Seurat 对象创建
# ------------------------------------------------------------------------------
read_loom_counts <- function(loom_file, sample, chunk_cells = 5000) {
  if (!file.exists(loom_file)) {
    stop("Loom file not found: ", loom_file)
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
    stop("Gene count does not match matrix columns in: ", loom_file)
  }
  if (length(barcodes) != n_cells) {
    stop("Cell barcode count does not match matrix rows in: ", loom_file)
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

sc_list <- list()

print("步骤1/4：读取 introns loom 矩阵并创建 Seurat 对象...")

for (sample in sample_ids) {
  loom_file <- file.path(data_dir, paste0(sample, ".loom"))

  tryCatch({
    counts <- read_loom_counts(loom_file, sample = sample)
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)

    if (ncol(sc_obj) == 0) next

    prefix <- substr(sample, 1, 1)
    suffix <- substr(sample, 2, 2)
    sc_obj$Tissue <- as.character(tissue_map[prefix])
    sc_obj$Group <- as.character(group_map[suffix])
    sc_obj$Condition <- ifelse(suffix == "2", "Cold", "NonCold")
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")

    sc_list[[sample]] <- sc_obj
    print(paste("   样本读取完成:", sample, "| 初始细胞数:", ncol(sc_obj)))

    rm(counts, sc_obj)
    gc(verbose = FALSE)
  }, error = function(e) {
    message(paste("   处理样本时出错", sample, ":", e$message))
  })
}

if (length(sc_list) == 0) {
  stop("没有样本成功读取，请检查 loom 输入路径。")
}

sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = names(sc_list))
sc_combined <- JoinLayers(sc_combined)

# ------------------------------------------------------------------------------
# 3. 综合质控：数值过滤 + 靶向污染过滤 + 核糖体基因过滤
# ------------------------------------------------------------------------------
print("步骤2/4：执行 QC 过滤和核糖体基因移除...")

sc_combined <- subset(sc_combined, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 15)

counts_mat <- GetAssayData(sc_combined, layer = "counts")
cells_keep <- colnames(sc_combined)

# 与当前参考脚本保持一致：血小板/红细胞靶向过滤逻辑保留在注释中，
# 但默认不启用。
# if ("Pf4" %in% rownames(counts_mat)) {
#   cells_keep <- cells_keep[counts_mat["Pf4", cells_keep] == 0]
# }
# for (rbc_gene in c("Hbb-bs", "Hba-a1")) {
#   if (rbc_gene %in% rownames(counts_mat)) {
#     cells_keep <- cells_keep[counts_mat[rbc_gene, cells_keep] <= 1]
#   }
# }

sc_combined <- subset(sc_combined, cells = cells_keep)
print(paste("   血小板/红细胞过滤步骤后细胞数:", ncol(sc_combined)))

sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = "^Rp[sl]")
ribo_genes <- grep("^Rp[sl]", rownames(sc_combined), value = TRUE, ignore.case = TRUE)
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]
sc_combined <- subset(sc_combined, features = non_ribo_genes)
print(paste("   移除核糖体基因后基因数:", nrow(sc_combined)))

# ------------------------------------------------------------------------------
# 4. 按组织拆分并分别运行 scDblFinder
# ------------------------------------------------------------------------------
print("步骤3/4：按组织拆分并运行 scDblFinder...")

sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name) {
  print(paste(">>> 正在处理组织:", tissue_name, "| QC 后细胞数:", ncol(obj)))
  if (ncol(obj) < 50) {
    return(NULL)
  }

  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5)

  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class

    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    print(paste("   [scDblFinder] 双细胞数:", n_dbl, "(", round(n_dbl / ncol(obj) * 100, 2), "%)"))

    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    print(paste("   [Filter] 去双胞后细胞数:", ncol(obj)))
  }, error = function(e) {
    message(paste("   scDblFinder 警告:", e$message, "- 保留所有细胞继续"))
  })

  obj
}

sc_by_tissue_clean <- lapply(names(sc_by_tissue), function(x) {
  run_standard_pipeline(sc_by_tissue[[x]], x)
})
names(sc_by_tissue_clean) <- names(sc_by_tissue)
sc_by_tissue_clean <- sc_by_tissue_clean[!sapply(sc_by_tissue_clean, is.null)]

if (length(sc_by_tissue_clean) == 0) {
  stop("QC/去双胞后没有剩余组织对象。")
}

# ------------------------------------------------------------------------------
# 5. 合并并保存
# ------------------------------------------------------------------------------
print("步骤4/4：合并清洗后的对象并保存 introns QC 输出...")

if (length(sc_by_tissue_clean) == 1) {
  sc_final <- sc_by_tissue_clean[[1]]
} else {
  sc_final <- merge(sc_by_tissue_clean[[1]], y = sc_by_tissue_clean[2:length(sc_by_tissue_clean)])
}
sc_final <- JoinLayers(sc_final)

qsave(sc_final, output_file)
print(paste("完成。最终高质量细胞数:", ncol(sc_final)))
print(paste("已保存:", output_file))
