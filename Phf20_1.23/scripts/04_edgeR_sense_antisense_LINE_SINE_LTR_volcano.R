#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20_1.23"
analysis_dir <- file.path(workdir, "sense_antisense_ERV")
results_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(results_dir, "metadata.tsv")
if (!file.exists(metadata_file)) stop("Metadata not found. Run 03_merge_TEcount_matrix.R first.")
metadata <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
metadata$group <- factor(metadata$group, levels = c("Scr", "Phf20"))

run_edger <- function(count_file, label) {
  if (!file.exists(count_file)) stop("Count file not found: ", count_file)
  dat <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("feature_id", "TE_class", metadata$sample) %in% colnames(dat))) {
    stop("Input file missing required columns: ", count_file)
  }

  count_mat <- as.matrix(dat[, metadata$sample, drop = FALSE])
  rownames(count_mat) <- make.unique(dat$feature_id)
  storage.mode(count_mat) <- "integer"

  y <- DGEList(counts = count_mat, group = metadata$group)
  keep <- filterByExpr(y, group = metadata$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ group, data = metadata)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupPhf20")

  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tested <- dat[keep, , drop = FALSE]
  idx <- match(rownames(res), rownames(y$counts))
  res$feature_id <- tested$feature_id[idx]
  res$TE_class <- tested$TE_class[idx]
  res$mean_count <- rowMeans(count_mat[keep, , drop = FALSE])[idx]
  res <- res[, c("feature_id", "TE_class", "logFC", "logCPM", "F", "PValue", "FDR", "mean_count")]

  res$regulation <- "Not significant"
  res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC >= 1] <- "Up in Phf20"
  res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC <= -1] <- "Down in Phf20"
  res$regulation <- factor(res$regulation, levels = c("Up in Phf20", "Down in Phf20", "Not significant"))

  cpm_mat <- cpm(y, log = TRUE, prior.count = 2)
  rownames(cpm_mat) <- res$feature_id[match(rownames(cpm_mat), rownames(y$counts))]

  write.table(res, file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_results.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(cbind(feature_id = rownames(cpm_mat), as.data.frame(cpm_mat, check.names = FALSE)),
              file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_logCPM_matrix.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  class_summary <- aggregate(
    cbind(
      tested_features = rep(1, nrow(res)),
      up_in_Phf20_FDR_0.05_logFC_1 = as.integer(res$regulation == "Up in Phf20"),
      down_in_Phf20_FDR_0.05_logFC_1 = as.integer(res$regulation == "Down in Phf20")
    ) ~ TE_class,
    data = res,
    FUN = sum
  )
  total <- data.frame(
    TE_class = "LINE_SINE_LTR",
    tested_features = nrow(res),
    up_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Up in Phf20"),
    down_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Down in Phf20")
  )
  summary_df <- rbind(total, class_summary)
  write.table(summary_df, file.path(results_dir, paste0("edgeR_", label, "_LINE_SINE_LTR_summary.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  plot_df <- res[is.finite(res$logFC) & !is.na(res$PValue), , drop = FALSE]
  plot_df$neg_log10_fdr <- -log10(pmax(plot_df$FDR, .Machine$double.xmin))
  set.seed(20260617)
  y_jitter_sd <- max(0.02, 0.006 * diff(range(plot_df$neg_log10_fdr, na.rm = TRUE)))
  plot_df$jitter_logFC <- plot_df$logFC + rnorm(nrow(plot_df), 0, 0.035)
  plot_df$jitter_neg_log10_fdr <- pmax(0, plot_df$neg_log10_fdr + rnorm(nrow(plot_df), 0, y_jitter_sd))
  plot_df$label <- ""
  sig_df <- plot_df[plot_df$regulation != "Not significant", , drop = FALSE]
  if (nrow(sig_df) > 0) {
    sig_idx <- order(sig_df$FDR, -abs(sig_df$logFC))[seq_len(min(12, nrow(sig_df)))]
    plot_df$label[match(sig_df$feature_id[sig_idx], plot_df$feature_id)] <- sig_df$feature_id[sig_idx]
  }

  count_label <- paste0(
    "Up in Phf20: ", sum(plot_df$regulation == "Up in Phf20"),
    "\nDown in Phf20: ", sum(plot_df$regulation == "Down in Phf20"),
    "\nTested: ", nrow(plot_df)
  )

  p <- ggplot(plot_df, aes(jitter_logFC, jitter_neg_log10_fdr, color = regulation, shape = TE_class)) +
    geom_point(alpha = 0.75, size = 1.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, min.segment.length = 0) +
    annotate("label",
             x = max(plot_df$logFC, na.rm = TRUE) - 0.04 * diff(range(plot_df$logFC, na.rm = TRUE)),
             y = max(plot_df$neg_log10_fdr, na.rm = TRUE) * 0.95,
             label = count_label, hjust = 1, vjust = 1, size = 3.4, fill = "white") +
    scale_color_manual(values = c("Up in Phf20" = "#D73027", "Down in Phf20" = "#4575B4", "Not significant" = "grey70"),
                       drop = FALSE) +
    scale_shape_manual(values = c("LINE" = 16, "SINE" = 17, "LTR" = 15)) +
    labs(title = paste0(label, " LINE/SINE/LTR differential expression"),
         subtitle = "edgeR QL GLM; logFC = Phf20 / Scr; significant: FDR < 0.05 and |logFC| >= 1",
         x = "log2 fold change", y = "-log10(FDR)", color = NULL, shape = "TE class") +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(face = "bold"), legend.position = "right")

  ggsave(file.path(plots_dir, paste0("volcano_", label, "_LINE_SINE_LTR.pdf")), p, width = 7.5, height = 6)
  ggsave(file.path(plots_dir, paste0("volcano_", label, "_LINE_SINE_LTR.png")), p, width = 7.5, height = 6, dpi = 300)

  list(res = res, summary = summary_df)
}

sense <- run_edger(file.path(results_dir, "TE_sense_LINE_SINE_LTR_counts.tsv"), "sense")
antisense <- run_edger(file.path(results_dir, "TE_antisense_LINE_SINE_LTR_counts.tsv"), "antisense")

message("Sense summary:")
print(sense$summary, row.names = FALSE)
message("Antisense summary:")
print(antisense$summary, row.names = FALSE)
message("Output directory: ", results_dir)
message("Plot directory: ", plots_dir)
