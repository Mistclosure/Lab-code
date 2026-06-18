#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
counts_dir <- file.path(workdir, "counts_antisense_TEcount")
outdir <- file.path(workdir, "antisense_TEcount_LINE_SINE_LTR_analysis")
files_dir <- file.path(outdir, "files")
plots_dir <- file.path(outdir, "plots")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

prefix <- "Phf20_ribo_deleption_antisense_TEcount"
target_classes <- c("LINE", "SINE", "LTR")

cnt_files <- list.files(
  counts_dir,
  pattern = "_antisense\\.cntTable$",
  full.names = TRUE
)

if (length(cnt_files) == 0) {
  stop("No antisense .cntTable files found in: ", counts_dir)
}

cnt_files <- sort(cnt_files)

sample_names <- sub("_antisense\\.cntTable$", "", basename(cnt_files))
sample_names <- make.names(sample_names, unique = TRUE)

read_cnt <- function(path, sample_name) {
  dat <- read.delim(
    path,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    quote = "\""
  )
  if (ncol(dat) != 2) {
    stop("Unexpected column number in ", path, ": ", ncol(dat))
  }
  colnames(dat) <- c("feature_id", sample_name)
  dat[[sample_name]] <- as.numeric(dat[[sample_name]])
  dat
}

count_list <- Map(read_cnt, cnt_files, sample_names)
merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), count_list)
merged[is.na(merged)] <- 0

count_cols <- setdiff(colnames(merged), "feature_id")
for (col in count_cols) {
  merged[[col]] <- as.numeric(merged[[col]])
}

feature_type <- ifelse(grepl("^ENSMUSG", merged$feature_id), "gene", "TE")
merged_with_type <- cbind(
  feature_id = merged$feature_id,
  feature_type = feature_type,
  merged[, count_cols, drop = FALSE]
)

write.csv(
  merged_with_type,
  file.path(files_dir, paste0(prefix, "_merged_gene_TE_counts.csv")),
  quote = FALSE,
  row.names = FALSE
)

te_df <- merged_with_type[merged_with_type$feature_type == "TE", , drop = FALSE]
te_df$TE_class <- sub("^.*:([^:]+)$", "\\1", te_df$feature_id)
target_df <- te_df[te_df$TE_class %in% target_classes, , drop = FALSE]

if (nrow(target_df) == 0) {
  stop("No LINE/SINE/LTR TE features found in antisense merged counts.")
}

target_counts_df <- target_df[, c("feature_id", "TE_class", count_cols), drop = FALSE]
write.csv(
  target_counts_df,
  file.path(files_dir, paste0(prefix, "_LINE_SINE_LTR_counts.csv")),
  quote = FALSE,
  row.names = FALSE
)

group <- ifelse(grepl("Phf20", count_cols, ignore.case = TRUE), "Phf20", "Scr")
group <- factor(group, levels = c("Scr", "Phf20"))

sample_metadata <- data.frame(
  sample = count_cols,
  group = group,
  stringsAsFactors = FALSE
)
write.csv(
  sample_metadata,
  file.path(files_dir, paste0(prefix, "_sample_metadata.csv")),
  quote = FALSE,
  row.names = FALSE
)

count_mat <- as.matrix(target_df[, count_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- make.unique(target_df$feature_id)

y <- DGEList(counts = count_mat, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "groupPhf20")

res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
matched_idx <- match(rownames(res), rownames(y$counts))
tested_target <- target_df[keep, , drop = FALSE]
res$feature_id <- tested_target$feature_id[matched_idx]
res$TE_class <- tested_target$TE_class[matched_idx]
res$mean_count <- rowMeans(count_mat[keep, , drop = FALSE])[matched_idx]
res <- res[, c("feature_id", "TE_class", "logFC", "logCPM", "F", "PValue", "FDR", "mean_count")]

res$regulation <- "Not significant"
res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC >= 1] <- "Up in Phf20"
res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC <= -1] <- "Down in Phf20"
res$regulation <- factor(
  res$regulation,
  levels = c("Up in Phf20", "Down in Phf20", "Not significant")
)

write.csv(
  res,
  file.path(files_dir, paste0(prefix, "_edgeR_DE_LINE_SINE_LTR.csv")),
  quote = FALSE,
  row.names = FALSE
)

class_summary <- data.frame(
  TE_class = target_classes,
  total_features = vapply(target_classes, function(cls) sum(target_df$TE_class == cls), numeric(1)),
  tested_features = vapply(target_classes, function(cls) sum(res$TE_class == cls), numeric(1)),
  up_in_Phf20_FDR_0.05_logFC_1 = vapply(
    target_classes,
    function(cls) sum(res$TE_class == cls & res$regulation == "Up in Phf20"),
    numeric(1)
  ),
  down_in_Phf20_FDR_0.05_logFC_1 = vapply(
    target_classes,
    function(cls) sum(res$TE_class == cls & res$regulation == "Down in Phf20"),
    numeric(1)
  ),
  stringsAsFactors = FALSE
)
summary_df <- rbind(
  data.frame(
    TE_class = "LINE_SINE_LTR",
    total_features = nrow(target_df),
    tested_features = nrow(res),
    up_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Up in Phf20"),
    down_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Down in Phf20"),
    stringsAsFactors = FALSE
  ),
  class_summary
)
write.csv(
  summary_df,
  file.path(files_dir, paste0(prefix, "_edgeR_DE_LINE_SINE_LTR_summary.csv")),
  quote = FALSE,
  row.names = FALSE
)

plot_df <- res[is.finite(res$logFC) & !is.na(res$PValue), , drop = FALSE]
plot_df$neg_log10_fdr <- -log10(pmax(plot_df$FDR, .Machine$double.xmin))

set.seed(20260617)
x_jitter_sd <- 0.035
y_range <- diff(range(plot_df$neg_log10_fdr, na.rm = TRUE))
y_jitter_sd <- max(0.02, 0.006 * y_range)
plot_df$jitter_logFC <- plot_df$logFC + rnorm(nrow(plot_df), mean = 0, sd = x_jitter_sd)
plot_df$jitter_neg_log10_fdr <- pmax(
  0,
  plot_df$neg_log10_fdr + rnorm(nrow(plot_df), mean = 0, sd = y_jitter_sd)
)

plot_df$label <- ""
sig_df <- plot_df[plot_df$regulation != "Not significant", , drop = FALSE]
if (nrow(sig_df) > 0) {
  sig_idx <- order(sig_df$FDR, -abs(sig_df$logFC))[seq_len(min(12, nrow(sig_df)))]
  labeled_features <- sig_df$feature_id[sig_idx]
  plot_df$label[match(labeled_features, plot_df$feature_id)] <- labeled_features
}

up_n <- sum(plot_df$regulation == "Up in Phf20")
down_n <- sum(plot_df$regulation == "Down in Phf20")
tested_n <- nrow(plot_df)

x_range <- range(plot_df$logFC, na.rm = TRUE)
x_pos <- max(plot_df$logFC, na.rm = TRUE) - 0.04 * diff(x_range)
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
    title = "Antisense LINE/SINE/LTR differential expression",
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

pdf_file <- file.path(plots_dir, paste0(prefix, "_edgeR_volcano_LINE_SINE_LTR.pdf"))
png_file <- file.path(plots_dir, paste0(prefix, "_edgeR_volcano_LINE_SINE_LTR.png"))
ggsave(pdf_file, p, width = 7.5, height = 6)
ggsave(png_file, p, width = 7.5, height = 6, dpi = 300)

message("Merged antisense files: ", length(cnt_files))
message("Total merged features: ", nrow(merged_with_type))
message("LINE/SINE/LTR features before filter: ", nrow(target_df))
message("LINE/SINE/LTR features tested: ", tested_n)
message("Up in Phf20: ", up_n)
message("Down in Phf20: ", down_n)
message("Output files directory: ", files_dir)
message("Output plots directory: ", plots_dir)
print(sample_metadata, row.names = FALSE)
print(summary_df, row.names = FALSE)
