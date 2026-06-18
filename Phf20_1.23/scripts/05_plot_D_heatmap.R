#!/usr/bin/env Rscript

if (!requireNamespace("pheatmap", quietly = TRUE)) {
  stop("pheatmap is not installed in the current R environment.")
}

suppressPackageStartupMessages({
  library(pheatmap)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20_1.23"
analysis_dir <- file.path(workdir, "sense_antisense_ERV")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

metadata <- read.delim(file.path(results_dir, "metadata.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
metadata$group <- factor(metadata$group, levels = c("Scr", "Phf20"))
ann_row <- data.frame(group = metadata$group)
rownames(ann_row) <- metadata$sample

plot_one <- function(label, title) {
  res_file <- file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_results.tsv"))
  mat_file <- file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_logCPM_matrix.tsv"))
  if (!file.exists(res_file) || !file.exists(mat_file)) {
    stop("Missing edgeR result or logCPM matrix for ", label, ". Run 04_edgeR_sense_antisense_LINE_SINE_LTR_volcano.R first.")
  }

  res <- read.delim(res_file, check.names = FALSE, stringsAsFactors = FALSE)
  mat_df <- read.delim(mat_file, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(mat_df) <- mat_df$feature_id
  mat <- as.matrix(mat_df[, metadata$sample, drop = FALSE])

  sig <- res[!is.na(res$FDR) & res$FDR < 0.05 & abs(res$logFC) > 1, , drop = FALSE]
  used_fallback <- FALSE
  if (nrow(sig) == 0) {
    used_fallback <- TRUE
    sig <- head(res[order(res$FDR, res$PValue), , drop = FALSE], 50)
    message("No significant ", label, " LINE/SINE/LTR by FDR < 0.05 and |logFC| > 1; using top 50 by FDR for heatmap.")
  }

  write.table(sig, file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_heatmap_features.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  keep_ids <- intersect(sig$feature_id, rownames(mat))
  if (length(keep_ids) == 0) stop("No heatmap features found in logCPM matrix for ", label)
  plot_mat <- mat[keep_ids, , drop = FALSE]
  z <- t(scale(t(plot_mat)))
  z[!is.finite(z)] <- 0
  z_t <- t(z)
  cluster_features <- ncol(z_t) >= 2

  heatmap_title <- title
  if (used_fallback) heatmap_title <- paste0(title, " (top 50 by FDR)")

  pdf_file <- file.path(plots_dir, paste0("plot_D_", label, "_LINE_SINE_LTR_heatmap.pdf"))
  png_file <- file.path(plots_dir, paste0("plot_D_", label, "_LINE_SINE_LTR_heatmap.png"))

  pdf(pdf_file, width = 10, height = 4.8)
  pheatmap(z_t, cluster_rows = FALSE, cluster_cols = cluster_features,
           annotation_row = ann_row, show_colnames = FALSE,
           main = heatmap_title, color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()

  png(png_file, width = 3000, height = 1440, res = 300)
  pheatmap(z_t, cluster_rows = FALSE, cluster_cols = cluster_features,
           annotation_row = ann_row, show_colnames = FALSE,
           main = heatmap_title, color = colorRampPalette(c("blue", "white", "red"))(100))
  dev.off()

  message("Wrote: ", pdf_file)
  message("Wrote: ", png_file)
}

plot_one("sense", "Differential expression of sense stranded LINE/SINE/LTR")
plot_one("antisense", "Differential expression of anti-sense stranded LINE/SINE/LTR")
