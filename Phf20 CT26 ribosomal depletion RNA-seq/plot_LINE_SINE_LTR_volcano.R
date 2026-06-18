#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
files_dir <- file.path(workdir, "files")
plots_dir <- file.path(workdir, "plots")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

prefix <- "Phf20_ribo_deleption"
de_file <- file.path(files_dir, paste0(prefix, "_edgeR_DE_TE.csv"))

if (!file.exists(de_file)) {
  stop("TE differential result file not found: ", de_file)
}

de <- read.csv(de_file, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
required_cols <- c("feature_id", "feature_type", "logFC", "PValue", "FDR", "regulation")
if (!all(required_cols %in% colnames(de))) {
  stop("Input file is missing required columns: ",
       paste(setdiff(required_cols, colnames(de)), collapse = ", "))
}

de$TE_class <- sub("^.*:([^:]+)$", "\\1", de$feature_id)
target_classes <- c("LINE", "SINE", "LTR")
plot_df <- de[de$TE_class %in% target_classes, , drop = FALSE]

if (nrow(plot_df) == 0) {
  stop("No LINE/SINE/LTR records found in: ", de_file)
}

plot_df <- plot_df[is.finite(plot_df$logFC) & !is.na(plot_df$PValue), , drop = FALSE]
plot_df$neg_log10_fdr <- -log10(pmax(plot_df$FDR, .Machine$double.xmin))
plot_df$regulation <- "Not significant"
plot_df$regulation[!is.na(plot_df$FDR) & plot_df$FDR < 0.05 & plot_df$logFC >= 1] <- "Up in Phf20"
plot_df$regulation[!is.na(plot_df$FDR) & plot_df$FDR < 0.05 & plot_df$logFC <= -1] <- "Down in Phf20"
plot_df$regulation <- factor(
  plot_df$regulation,
  levels = c("Up in Phf20", "Down in Phf20", "Not significant")
)

# Gaussian jitter is used only for plotting to reduce overplotting.
# Differential statistics written to CSV remain unchanged.
set.seed(20260617)
x_jitter_sd <- 0.035
y_jitter_sd <- max(0.02, 0.006 * diff(range(plot_df$neg_log10_fdr, na.rm = TRUE)))
plot_df$jitter_logFC <- plot_df$logFC + rnorm(nrow(plot_df), mean = 0, sd = x_jitter_sd)
plot_df$jitter_neg_log10_fdr <- pmax(
  0,
  plot_df$neg_log10_fdr + rnorm(nrow(plot_df), mean = 0, sd = y_jitter_sd)
)

up_n <- sum(plot_df$regulation == "Up in Phf20")
down_n <- sum(plot_df$regulation == "Down in Phf20")
tested_n <- nrow(plot_df)

class_summary <- as.data.frame(table(plot_df$TE_class), stringsAsFactors = FALSE)
colnames(class_summary) <- c("TE_class", "tested_features")
class_summary$up_in_Phf20_FDR_0.05_logFC_1 <- vapply(
  class_summary$TE_class,
  function(cls) sum(plot_df$TE_class == cls & plot_df$regulation == "Up in Phf20"),
  numeric(1)
)
class_summary$down_in_Phf20_FDR_0.05_logFC_1 <- vapply(
  class_summary$TE_class,
  function(cls) sum(plot_df$TE_class == cls & plot_df$regulation == "Down in Phf20"),
  numeric(1)
)

write.csv(
  plot_df,
  file.path(files_dir, paste0(prefix, "_edgeR_DE_TE_LINE_SINE_LTR.csv")),
  quote = FALSE,
  row.names = FALSE
)
write.csv(
  class_summary,
  file.path(files_dir, paste0(prefix, "_edgeR_DE_TE_LINE_SINE_LTR_summary.csv")),
  quote = FALSE,
  row.names = FALSE
)

plot_df$label <- ""
sig_df <- plot_df[plot_df$regulation != "Not significant", , drop = FALSE]
if (nrow(sig_df) > 0) {
  sig_idx <- order(sig_df$FDR, -abs(sig_df$logFC))[seq_len(min(12, nrow(sig_df)))]
  labeled_features <- sig_df$feature_id[sig_idx]
  plot_df$label[match(labeled_features, plot_df$feature_id)] <- labeled_features
}

x_pos <- max(plot_df$logFC, na.rm = TRUE) - 0.04 * diff(range(plot_df$logFC, na.rm = TRUE))
y_pos <- max(plot_df$neg_log10_fdr, na.rm = TRUE) * 0.95
count_label <- paste0(
  "Up in Phf20: ", up_n,
  "\nDown in Phf20: ", down_n,
  "\nTested LINE/SINE/LTR: ", tested_n
)

p <- ggplot(plot_df, aes(x = jitter_logFC, y = jitter_neg_log10_fdr, color = regulation, shape = TE_class)) +
  geom_point(alpha = 0.75, size = 1.5) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey40") +
  geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, min.segment.length = 0) +
  annotate(
    "label",
    x = x_pos,
    y = y_pos,
    label = count_label,
    hjust = 1,
    vjust = 1,
    size = 3.6,
    fill = "white"
  ) +
  scale_color_manual(
    values = c(
      "Up in Phf20" = "#D73027",
      "Down in Phf20" = "#4575B4",
      "Not significant" = "grey70"
    ),
    drop = FALSE
  ) +
  scale_shape_manual(values = c("LINE" = 16, "SINE" = 17, "LTR" = 15)) +
  labs(
    title = "LINE/SINE/LTR differential expression",
    subtitle = "edgeR QL GLM; logFC = Phf20 / Scr; significant: FDR < 0.05 and |logFC| >= 1",
    x = "log2 fold change",
    y = "-log10(FDR)",
    color = NULL,
    shape = "TE class"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

pdf_file <- file.path(plots_dir, paste0(prefix, "_edgeR_volcano_TE_LINE_SINE_LTR.pdf"))
png_file <- file.path(plots_dir, paste0(prefix, "_edgeR_volcano_TE_LINE_SINE_LTR.png"))
ggsave(pdf_file, p, width = 7.5, height = 6)
ggsave(png_file, p, width = 7.5, height = 6, dpi = 300)

message("Input TE features: ", nrow(de))
message("LINE/SINE/LTR features tested: ", tested_n)
message("Up in Phf20: ", up_n)
message("Down in Phf20: ", down_n)
message("Output plot PDF: ", pdf_file)
message("Output plot PNG: ", png_file)
print(class_summary, row.names = FALSE)
