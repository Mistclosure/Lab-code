# ==============================================================================
# 08_combined_macrophage_monocyte_monocle2_velocity_introns.R
# coldmouse_introns: Aorta Macrophage + PBMC Monocyte combined monocle2 + velocity exports
# ============================================================================== 
# 输入：qs/coldmouse_introns_sc_by_tissue_no_harmony_<cluster_method>.qs
# 输出：
#   1) 合并后的 Aorta Macrophage + PBMC Monocyte Seurat 对象
#   2) combined UMAP 坐标、metadata、loom obs_name 映射，用于 scVelo
#   3) 合并子集的一次 monocle2 拟时序结果
# 注意：RNA velocity 由 09_scvelo_velocity_combined_introns.py 使用 scVelo 执行。
# ============================================================================== 

PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))

config <- load_config()
setup_performance(config)

suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
  library(ggplot2)
  library(dplyr)
  library(Matrix)
  library(Biobase)
  library(monocle)
})

# ------------------------------------------------------------------------------
# monocle2 compatibility patches for current dplyr / igraph
# ------------------------------------------------------------------------------
for (fn_name in c("group_by_", "summarise_", "mutate_", "filter_", "arrange_",
                  "select_", "slice_", "distinct_")) {
  if (!exists(fn_name, envir = globalenv())) {
    base_fn <- gsub("_$", "", fn_name)
    fn_def <- eval(parse(text = paste0("function(.data, ...) { dplyr::", base_fn, "(.data, ...) }")))
    assign(fn_name, fn_def, envir = globalenv())
  }
}

if (requireNamespace("igraph", quietly = TRUE)) {
  orig_dfs <- igraph::dfs
  patched_dfs <- function(graph, root, neimode = NULL, mode = NULL, ...) {
    if (!is.null(neimode) && is.null(mode)) mode <- neimode
    orig_dfs(graph, root, mode = mode, ...)
  }
  environment(patched_dfs) <- environment(orig_dfs)
  assignInNamespace("dfs", patched_dfs, ns = "igraph")
}

patch_monocle_nei <- function() {
  ns <- getNamespace("monocle")
  for (fn_name in ls(ns)) {
    fn <- get(fn_name, envir = ns)
    if (!is.function(fn)) next
    body_txt <- paste(deparse(body(fn)), collapse = "\n")
    if (!grepl("nei(", body_txt, fixed = TRUE)) next
    patched_txt <- gsub("(?<![.])nei\\(", ".nei(", body_txt, perl = TRUE)
    patched_body <- tryCatch(parse(text = patched_txt)[[1]], error = function(e) NULL)
    if (is.null(patched_body)) next
    body(fn) <- patched_body
    tryCatch(assignInNamespace(fn_name, fn, ns = "monocle"), error = function(e) NULL)
  }
}
patch_monocle_nei()

make_compat_select_ <- function() {
  function(.data, ..., .dots = NULL) {
    lazies <- lazyeval::lazy_dots(...)
    vals <- lapply(lazies, function(z) {
      tryCatch(lazyeval::lazy_eval(z), error = function(e) as.character(z$expr))
    })
    if (!is.null(.dots) && length(.dots) > 0) vals <- c(vals, .dots)
    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0) return(.data)
    cols <- unlist(lapply(vals, function(x) {
      if (length(x) == 0) return(character(0))
      x1 <- x[[1]]
      if (is.numeric(x1)) return(names(.data)[as.integer(x1)])
      as.character(x1)
    }), use.names = FALSE)
    cols <- cols[nzchar(cols)]
    out <- dplyr::select(.data, dplyr::all_of(cols))
    new_names <- names(vals)
    if (!is.null(new_names)) {
      named <- !is.na(new_names) & nzchar(new_names)
      if (any(named)) {
        idx <- match(cols[named], names(out))
        ok <- !is.na(idx)
        names(out)[idx[ok]] <- new_names[named][ok]
      }
    }
    out
  }
}
make_compat_group_by_ <- function() {
  function(.data, ..., .dots = NULL) {
    lazies <- lazyeval::lazy_dots(...)
    vals <- lapply(lazies, function(z) {
      tryCatch(lazyeval::lazy_eval(z), error = function(e) as.character(z$expr))
    })
    if (!is.null(.dots) && length(.dots) > 0) vals <- c(vals, .dots)
    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0) return(dplyr::group_by(.data))
    cols <- unlist(lapply(vals, function(x) as.character(x[[1]])), use.names = FALSE)
    cols <- cols[nzchar(cols)]
    dplyr::group_by(.data, dplyr::across(dplyr::all_of(cols)))
  }
}
patch_binding <- function(env, name, value) {
  if (!exists(name, envir = env, inherits = FALSE)) return(invisible(FALSE))
  was_locked <- bindingIsLocked(name, env)
  if (was_locked) unlockBinding(name, env)
  assign(name, value, envir = env)
  if (was_locked) lockBinding(name, env)
  invisible(TRUE)
}
dplyr_compat_fns <- list(
  select_ = make_compat_select_(),
  group_by_ = make_compat_group_by_(),
  filter_ = function(.data, ..., .dots = NULL) dplyr::filter(.data, ...),
  arrange_ = function(.data, ..., .dots = NULL) dplyr::arrange(.data, ...),
  mutate_ = function(.data, ..., .dots = NULL) dplyr::mutate(.data, ...),
  summarise_ = function(.data, ..., .dots = NULL) dplyr::summarise(.data, ...),
  summarize_ = function(.data, ..., .dots = NULL) dplyr::summarise(.data, ...)
)
dplyr_ns <- getNamespace("dplyr")
for (fn_name in names(dplyr_compat_fns)) {
  patch_binding(dplyr_ns, fn_name, dplyr_compat_fns[[fn_name]])
  if ("package:dplyr" %in% search()) patch_binding(as.environment("package:dplyr"), fn_name, dplyr_compat_fns[[fn_name]])
  patch_binding(parent.env(getNamespace("monocle")), fn_name, dplyr_compat_fns[[fn_name]])
  assign(fn_name, dplyr_compat_fns[[fn_name]], envir = globalenv())
}

# ------------------------------------------------------------------------------
# Paths and parameters
# ------------------------------------------------------------------------------
output_tag <- "AortaMacrophage_PBMCMonocyte"
objects_dir <- get_results_dir(config, "objects")
rds_dir <- get_results_dir(config, "rds")
files_base_dir <- get_results_dir(config, "files")
plots_base_dir <- get_results_dir(config, "plots")
files_dir <- file.path(files_base_dir, "Combined", output_tag)
plots_dir <- file.path(plots_base_dir, "Combined", output_tag)
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

MAX_MONOCLE2_CELLS <- as.integer(Sys.getenv("COMBINED_MONOCLE2_MAX_CELLS", unset = "5000"))
MAX_ORDERING_GENES <- as.integer(Sys.getenv("COMBINED_MONOCLE2_MAX_ORDERING_GENES", unset = "1200"))
set.seed(20260604)

make_loom_obs_name <- function(cell_names) {
  sample_id <- sub("_.*$", "", cell_names)
  barcode <- sub("^[^_]+_", "", cell_names)
  barcode <- gsub("-", "_", barcode, fixed = TRUE)
  paste0(sample_id, ":", barcode)
}

stratified_downsample_cells <- function(obj, max_cells = MAX_MONOCLE2_CELLS) {
  if (ncol(obj) <= max_cells) return(colnames(obj))
  meta <- obj@meta.data
  meta$.cell <- rownames(meta)
  meta$.stratum <- paste(meta$source_subset, meta$Group, meta$seurat_clusters, sep = "__")
  strata <- split(meta$.cell, meta$.stratum)
  strata <- strata[lengths(strata) > 0]
  target <- pmax(1, floor(lengths(strata) / ncol(obj) * max_cells))
  overflow <- sum(target) - max_cells
  if (overflow > 0) {
    reducible <- which(target > 1)
    while (overflow > 0 && length(reducible) > 0) {
      i <- reducible[which.max(target[reducible])]
      target[i] <- target[i] - 1
      overflow <- overflow - 1
      reducible <- which(target > 1)
    }
  }
  sampled <- unlist(Map(function(cells, n) sample(cells, min(length(cells), n)), strata, target), use.names = FALSE)
  if (length(sampled) > max_cells) sampled <- sample(sampled, max_cells)
  sampled
}

# ------------------------------------------------------------------------------
# Read and subset Seurat objects
# ------------------------------------------------------------------------------
sc_by_tissue_file <- file.path(objects_dir, paste0("coldmouse_introns_sc_by_tissue_no_harmony_", config$clustering$method, ".qs"))
cat("读取:", sc_by_tissue_file, "\n")
sc_by_tissue <- qread(sc_by_tissue_file)

stopifnot("Aorta" %in% names(sc_by_tissue), "PBMC" %in% names(sc_by_tissue))

aorta <- sc_by_tissue[["Aorta"]]
pbmc <- sc_by_tissue[["PBMC"]]
if (!"SingleR.labels" %in% colnames(aorta@meta.data) || !"SingleR.labels" %in% colnames(pbmc@meta.data)) {
  stop("SingleR.labels 字段不存在，请先完成 SingleR 注释。")
}

aorta_mac <- subset(aorta, subset = SingleR.labels == "Macrophages")
pbmc_mono <- subset(pbmc, subset = SingleR.labels == "Monocytes")
cat("Aorta Macrophages:", ncol(aorta_mac), "cells\n")
cat("PBMC Monocytes:", ncol(pbmc_mono), "cells\n")
if (ncol(aorta_mac) == 0 || ncol(pbmc_mono) == 0) stop("子集为空，检查 SingleR.labels。")

aorta_mac$source_subset <- "Aorta_Macrophage"
pbmc_mono$source_subset <- "PBMC_Monocyte"
aorta_mac$source_tissue <- "Aorta"
pbmc_mono$source_tissue <- "PBMC"
aorta_mac$source_celltype <- "Macrophage"
pbmc_mono$source_celltype <- "Monocyte"
aorta_mac$original_seurat_clusters <- as.character(aorta_mac$seurat_clusters)
pbmc_mono$original_seurat_clusters <- as.character(pbmc_mono$seurat_clusters)

aorta_mac <- JoinLayers(aorta_mac)
pbmc_mono <- JoinLayers(pbmc_mono)
combined <- merge(aorta_mac, y = pbmc_mono, project = "coldmouse_introns_AortaMacrophage_PBMCMonocyte")
combined$sample_id <- sub("_.*$", "", colnames(combined))
combined$loom_obs_name <- make_loom_obs_name(colnames(combined))
combined$seurat_cell <- colnames(combined)

cat("合并后细胞数:", ncol(combined), " genes:", nrow(combined), "\n")
print(table(combined$source_subset, combined$Group))

# ------------------------------------------------------------------------------
# Recompute combined UMAP for the merged subset
# ------------------------------------------------------------------------------
DefaultAssay(combined) <- "RNA"
combined <- NormalizeData(combined, verbose = FALSE)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = config$clustering$nfeatures, verbose = FALSE)
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = min(50, config$clustering$n_pcs), verbose = FALSE)
use_dims <- seq_len(min(30, ncol(Embeddings(combined, "pca"))))
combined <- FindNeighbors(combined, dims = use_dims, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, algorithm = 4, verbose = FALSE)
combined <- RunUMAP(combined, dims = use_dims, verbose = FALSE)
combined$combined_clusters <- as.character(combined$seurat_clusters)

combined_qs <- file.path(objects_dir, paste0("introns_combined_", output_tag, "_seurat.qs"))
qsave(combined, combined_qs)
cat("保存 combined Seurat:", combined_qs, "\n")

umap <- Embeddings(combined, "umap")
meta <- combined@meta.data
velocity_meta <- data.frame(
  seurat_cell = colnames(combined),
  loom_obs_name = meta$loom_obs_name,
  sample_id = meta$sample_id,
  source_subset = meta$source_subset,
  source_tissue = meta$source_tissue,
  source_celltype = meta$source_celltype,
  Tissue = meta$Tissue,
  Group = meta$Group,
  Condition = meta$Condition,
  SingleR.labels = meta$SingleR.labels,
  original_seurat_clusters = meta$original_seurat_clusters,
  combined_clusters = meta$combined_clusters,
  nCount_RNA = meta$nCount_RNA,
  nFeature_RNA = meta$nFeature_RNA,
  percent.mt = meta$percent.mt,
  umap_1 = umap[, 1],
  umap_2 = umap[, 2],
  stringsAsFactors = FALSE
)
metadata_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_velocity_metadata.csv"))
write.csv(velocity_meta, metadata_file, row.names = FALSE)
cat("保存 velocity metadata:", metadata_file, "\n")

sample_map <- data.frame(sample_id = sort(unique(velocity_meta$sample_id)), stringsAsFactors = FALSE)
sample_map <- data.frame(
  sample_id = sample_map$sample_id,
  loom_file = file.path(config$paths$rawdata_dir, paste0(sample_map$sample_id, ".loom")),
  stringsAsFactors = FALSE
)
loom_map_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_loom_file_map.csv"))
write.csv(sample_map, loom_map_file, row.names = FALSE)
cat("保存 loom map:", loom_map_file, "\n")

plot_umap <- function(color_by, filename, width = 8, height = 6) {
  p <- DimPlot(combined, reduction = "umap", group.by = color_by, pt.size = 0.25) +
    ggtitle(paste("Combined Aorta Macrophage + PBMC Monocyte:", color_by)) +
    theme_classic(base_size = 12)
  ggsave(file.path(plots_dir, paste0(filename, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(plots_dir, paste0(filename, ".pdf")), p, width = width, height = height)
}
plot_umap("source_subset", paste0("introns_combined_", output_tag, "_UMAP_source_subset"))
plot_umap("Group", paste0("introns_combined_", output_tag, "_UMAP_Group"))
plot_umap("combined_clusters", paste0("introns_combined_", output_tag, "_UMAP_combined_clusters"), width = 9, height = 6)

# ------------------------------------------------------------------------------
# Combined monocle2 trajectory
# ------------------------------------------------------------------------------
monocle_cells <- stratified_downsample_cells(combined, MAX_MONOCLE2_CELLS)
if (length(monocle_cells) < ncol(combined)) {
  cat("monocle2 分层抽样:", ncol(combined), "->", length(monocle_cells), "cells\n")
}
mono_obj <- subset(combined, cells = monocle_cells)
mono_obj <- JoinLayers(mono_obj)

counts <- GetAssayData(mono_obj, assay = "RNA", layer = "counts")
pd <- new("AnnotatedDataFrame", data = mono_obj@meta.data)
fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts), stringsAsFactors = FALSE)
fd <- new("AnnotatedDataFrame", data = fd)
cds <- newCellDataSet(
  as(counts, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  lowerDetectionLimit = 0.5,
  expressionFamily = negbinomial.size()
)

cat("estimateSizeFactors...\n")
cds <- tryCatch(estimateSizeFactors(cds), error = function(e) {
  cat("estimateSizeFactors 失败，手动计算:", e$message, "\n")
  cell_totals <- Matrix::colSums(exprs(cds))
  pData(cds)$Size_Factor <- cell_totals / median(cell_totals)
  cds
})
cat("estimateDispersions...\n")
cds <- tryCatch(estimateDispersions(cds), error = function(e) {
  cat("estimateDispersions 失败，继续:", e$message, "\n")
  cds
})

Idents(mono_obj) <- "combined_clusters"
all_markers <- FindAllMarkers(mono_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
markers_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_monocle2_cluster_markers.csv"))
write.csv(all_markers, markers_file, row.names = FALSE)

sig_markers <- all_markers %>% filter(p_val_adj < 0.05)
if (nrow(sig_markers) >= 100) {
  ordering_genes <- sig_markers %>% arrange(p_val_adj, desc(avg_log2FC)) %>% pull(gene) %>% unique()
} else {
  ordering_genes <- all_markers %>% group_by(cluster) %>% slice_max(n = 50, order_by = avg_log2FC, with_ties = FALSE) %>% arrange(desc(avg_log2FC)) %>% pull(gene) %>% unique()
}
ordering_genes <- ordering_genes[ordering_genes %in% rownames(cds)]
if (length(ordering_genes) > MAX_ORDERING_GENES) ordering_genes <- ordering_genes[seq_len(MAX_ORDERING_GENES)]
if (length(ordering_genes) < 100) stop("ordering genes < 100，无法稳定运行 monocle2。")
ordering_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_monocle2_ordering_genes.csv"))
write.csv(data.frame(gene = ordering_genes), ordering_file, row.names = FALSE)
cat("ordering genes:", length(ordering_genes), "\n")

cds <- setOrderingFilter(cds, ordering_genes)
cat("reduceDimension DDRTree...\n")
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cat("orderCells...\n")
cds <- orderCells(cds)

pseudotime_meta <- data.frame(
  seurat_cell = rownames(pData(cds)),
  loom_obs_name = pData(cds)$loom_obs_name,
  sample_id = pData(cds)$sample_id,
  source_subset = pData(cds)$source_subset,
  source_tissue = pData(cds)$source_tissue,
  source_celltype = pData(cds)$source_celltype,
  Group = pData(cds)$Group,
  combined_clusters = pData(cds)$combined_clusters,
  State = pData(cds)$State,
  Pseudotime = pData(cds)$Pseudotime,
  stringsAsFactors = FALSE
)
pseudotime_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_monocle2_pseudotime_metadata.csv"))
write.csv(pseudotime_meta, pseudotime_file, row.names = FALSE)

state_table <- as.data.frame.matrix(table(State = pData(cds)$State, subset = pData(cds)$source_subset))
state_table_file <- file.path(files_dir, paste0("introns_combined_", output_tag, "_monocle2_state_by_subset.csv"))
write.csv(state_table, state_table_file)

cds_qs <- file.path(objects_dir, paste0("introns_combined_", output_tag, "_monocle2_cds.qs"))
cds_rds <- file.path(rds_dir, paste0("introns_combined_", output_tag, "_monocle2_cds.rds"))
qsave(cds, cds_qs)
saveRDS(cds, cds_rds)

plot_configs <- list(
  list(color_by = "Pseudotime", w = 8, h = 6),
  list(color_by = "State", w = 8, h = 6),
  list(color_by = "source_subset", w = 8, h = 6),
  list(color_by = "Group", w = 8, h = 6),
  list(color_by = "combined_clusters", w = 10, h = 6)
)
for (pc in plot_configs) {
  p <- plot_cell_trajectory(cds, color_by = pc$color_by) +
    ggtitle(paste0("Combined Aorta Macrophage + PBMC Monocyte (", pc$color_by, ")"))
  ggsave(file.path(plots_dir, paste0("introns_combined_", output_tag, "_monocle2_trajectory_", pc$color_by, ".png")), p, width = pc$w, height = pc$h, dpi = 300)
  ggsave(file.path(plots_dir, paste0("introns_combined_", output_tag, "_monocle2_trajectory_", pc$color_by, ".pdf")), p, width = pc$w, height = pc$h)
}

cat("\n=== Combined monocle2 + velocity export 完成 ===\n")
cat("Seurat QS:", combined_qs, "\n")
cat("CDS QS:", cds_qs, "\n")
cat("Velocity metadata:", metadata_file, "\n")
cat("下一步 scVelo 脚本: 09_scvelo_velocity_combined_introns.py\n")
