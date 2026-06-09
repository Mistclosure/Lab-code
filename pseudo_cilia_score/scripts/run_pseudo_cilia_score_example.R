#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(qs)
})

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0) {
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
} else {
  getwd()
}
project_dir <- normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
if (!file.exists(file.path(project_dir, "R", "pseudo_cilia_score.R"))) {
  project_dir <- getwd()
}
source(file.path(project_dir, "R", "pseudo_cilia_score.R"))
setwd(project_dir)

# Edit these paths if your project stores inputs elsewhere.
seurat_path <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/qs/Seurat/Malignant_RNA_assay.qs"
cilia_gene_file <- "/mnt/disk1/qiuzerui/downloads/CRC/signature/ciliopathy_genes.csv"
if (!file.exists(cilia_gene_file)) {
  cilia_gene_file <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table2_primary_cilia_single_gene_confirmed.csv"
}

outdir <- "results/pseudo_cilia_score"
plot_dir <- file.path(outdir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(seurat_path)) {
  stop("Seurat object not found: ", seurat_path,
       "\nPlease edit seurat_path at the top of this script.")
}
if (!file.exists(cilia_gene_file)) {
  stop("Cilia gene file not found: ", cilia_gene_file,
       "\nPlease edit cilia_gene_file at the top of this script.")
}

obj <- qread(seurat_path)
cilia_df <- read.csv(cilia_gene_file, header = TRUE, check.names = FALSE)
gene_col <- intersect(c("Gene", "gene", "genes", "Symbol", "symbol"), colnames(cilia_df))[1]
if (is.na(gene_col)) {
  gene_col <- colnames(cilia_df)[1]
}
cilia_genes <- unique(as.character(cilia_df[[gene_col]]))

# Adapt these metadata columns to your Seurat object.
sample_col <- "orig.ident"
group_col <- "metastasis"
celltype_col <- NULL
subcluster_col <- "malignant_subcluster"
target_celltype <- NULL

# Example fallback for GSE132465: derive a group column if it is absent.
if (!group_col %in% colnames(obj@meta.data) && "TNM stage" %in% colnames(obj@meta.data)) {
  obj$metastasis <- ifelse(grepl("M0$", obj@meta.data[["TNM stage"]]), "Primary", "Metastasis")
}
if (!subcluster_col %in% colnames(obj@meta.data)) {
  subcluster_col <- NULL
}

result <- cal_pseudoCiliaScore(
  object = obj,
  cilia_genes = cilia_genes,
  assay = "RNA",
  layer = "counts",
  sample_col = sample_col,
  group_col = group_col,
  celltype_col = celltype_col,
  subcluster_col = subcluster_col,
  target_celltype = target_celltype,
  pseudo_mode = "chunk",
  cells_per_pseudo = 80,
  min_cells_per_pseudo = 30,
  min_total_count = 10,
  min_detect_pseudocells = 3,
  min_detect_rate = 0.01,
  remove_mt = TRUE,
  remove_ribo = FALSE,
  n.mean.bins = 20,
  n.detect.bins = 10,
  n.control = 50,
  num.rounds = 1000,
  seed = 123,
  return_matrices = FALSE,
  verbose = TRUE
)

save_pseudoCiliaScore_results(result, outdir)

score <- result$score
sample_summary <- result$sample_summary
gene_info <- result$gene_info
used_gene_info <- gene_info[gene_info$used_as_cilia_gene, , drop = FALSE]

p1 <- ggplot(score, aes(x = group, y = CiliaScore, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 1.2, alpha = 0.45) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(title = "Pseudo-cell CiliaScore", x = NULL, y = "CiliaScore")
ggsave(file.path(plot_dir, "pseudo_cell_CiliaScore_boxplot.png"), p1, width = 5.5, height = 4.5, dpi = 300)

p2 <- ggplot(sample_summary, aes(x = group, y = mean_CiliaScore, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 2.2) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  labs(title = "Sample-level mean CiliaScore", x = NULL, y = "Mean CiliaScore")
ggsave(file.path(plot_dir, "sample_mean_CiliaScore_boxplot.png"), p2, width = 5.5, height = 4.5, dpi = 300)

p3 <- ggplot(score, aes(x = background_score, y = target_score, color = group)) +
  geom_point(size = 1.6, alpha = 0.75) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Target score vs background score",
       x = "Background score", y = "Target score")
ggsave(file.path(plot_dir, "target_vs_background_score.png"), p3, width = 5.5, height = 4.5, dpi = 300)

p4 <- ggplot(used_gene_info, aes(x = detection_rate)) +
  geom_histogram(bins = 20, fill = "#4C78A8", color = "white") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(title = "Detection rate of used cilia genes",
       x = "Detection rate in pseudo-cells", y = "Gene count")
ggsave(file.path(plot_dir, "used_cilia_gene_detection_rate.png"), p4, width = 5.5, height = 4.5, dpi = 300)

if (length(unique(sample_summary$group)) == 2 && nrow(sample_summary) >= 4) {
  wt <- wilcox.test(mean_CiliaScore ~ group, data = sample_summary)
  writeLines(c(
    "Sample-level Wilcoxon rank-sum test:",
    paste0("p.value = ", signif(wt$p.value, 4)),
    "",
    "Interpretation note:",
    "Pseudo-cell plots are useful for visualization, but formal group comparison should use sample-level summaries."
  ), file.path(outdir, "sample_level_group_test.txt"))
}

message("Finished. Results saved to: ", normalizePath(outdir, mustWork = FALSE))
