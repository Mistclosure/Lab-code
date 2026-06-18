#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
analysis_dir <- file.path(workdir, "sense_antisense_ERV")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

sense_file <- file.path(results_dir, "edgeR_sense_LINE_SINE_LTR_results.tsv")
antisense_file <- file.path(results_dir, "edgeR_antisense_LINE_SINE_LTR_results.tsv")
if (!file.exists(sense_file) || !file.exists(antisense_file)) {
  stop("Missing edgeR result files. Run 04_edgeR_sense_antisense_LINE_SINE_LTR_volcano.R first.")
}

sense <- read.delim(sense_file, check.names = FALSE, stringsAsFactors = FALSE)
antisense <- read.delim(antisense_file, check.names = FALSE, stringsAsFactors = FALSE)

colnames(sense)[colnames(sense) != "feature_id"] <- paste0("sense_", colnames(sense)[colnames(sense) != "feature_id"])
colnames(antisense)[colnames(antisense) != "feature_id"] <- paste0("antisense_", colnames(antisense)[colnames(antisense) != "feature_id"])
merged <- merge(sense, antisense, by = "feature_id", all = FALSE)

merged$bidirectional_up <- !is.na(merged$sense_FDR) & !is.na(merged$antisense_FDR) &
  merged$sense_FDR < 0.05 & merged$antisense_FDR < 0.05 &
  merged$sense_logFC > 1 & merged$antisense_logFC > 1
merged$concordant_positive <- merged$sense_logFC > 0 & merged$antisense_logFC > 0
merged$mean_log2FC <- rowMeans(cbind(merged$sense_logFC, merged$antisense_logFC), na.rm = TRUE)
merged$fold_change_Phf20_vs_Scr <- 2 ^ merged$mean_log2FC
merged$combined_rank_score <- -log10(pmax(merged$sense_FDR, .Machine$double.xmin)) +
  -log10(pmax(merged$antisense_FDR, .Machine$double.xmin)) +
  abs(merged$mean_log2FC)

write.table(merged, file.path(results_dir, "edgeR_sense_antisense_LINE_SINE_LTR_merged_results.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_df <- merged[merged$bidirectional_up, , drop = FALSE]
subtitle <- "Significant bidirectional up LINE/SINE/LTR: sense and antisense FDR < 0.05, logFC > 1"
if (nrow(plot_df) == 0) {
  plot_df <- merged[merged$concordant_positive, , drop = FALSE]
  plot_df <- plot_df[order(-plot_df$combined_rank_score), , drop = FALSE]
  subtitle <- "No significant bidirectional up features; showing top concordant positive LINE/SINE/LTR"
}
plot_df <- head(plot_df[order(-plot_df$fold_change_Phf20_vs_Scr), , drop = FALSE], 10)
if (nrow(plot_df) == 0) stop("No concordant positive sense/antisense LINE/SINE/LTR available for plot E summary.")

plot_df$feature_id <- factor(plot_df$feature_id, levels = plot_df$feature_id)

p <- ggplot(plot_df, aes(x = feature_id, y = fold_change_Phf20_vs_Scr)) +
  geom_col(fill = "black", width = 0.72) +
  labs(
    title = "Bi-directional LINE/SINE/LTR transcripts",
    subtitle = subtitle,
    x = NULL,
    y = "Fold change (Phf20/Scr)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(size = 10, angle = 50, hjust = 1, vjust = 1)
  )

pdf_file <- file.path(plots_dir, "plot_E_bidirectional_LINE_SINE_LTR_fold_change.pdf")
png_file <- file.path(plots_dir, "plot_E_bidirectional_LINE_SINE_LTR_fold_change.png")
ggsave(pdf_file, p, width = 9.5, height = 4.8)
ggsave(png_file, p, width = 9.5, height = 4.8, dpi = 300)

message("Wrote: ", pdf_file)
message("Wrote: ", png_file)
message("Plotted features: ", nrow(plot_df))
