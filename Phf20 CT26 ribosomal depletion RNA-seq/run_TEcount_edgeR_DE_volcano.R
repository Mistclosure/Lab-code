#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
files_dir <- file.path(workdir, "files")
plots_dir <- file.path(workdir, "plots")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

prefix <- "Phf20_ribo_deleption"
count_file <- file.path(files_dir, paste0(prefix, "_TEcount_merged_gene_TE_counts.csv"))

if (!file.exists(count_file)) {
  stop("Count file not found: ", count_file)
}

counts_df <- read.csv(
  count_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_cols <- c("feature_id", "feature_type")
if (!all(required_cols %in% colnames(counts_df))) {
  stop("Input file must contain columns: ", paste(required_cols, collapse = ", "))
}

sample_cols <- setdiff(colnames(counts_df), required_cols)
if (length(sample_cols) < 4) {
  stop("Expected at least 4 sample columns, found: ", length(sample_cols))
}

group <- ifelse(grepl("Phf20", sample_cols, ignore.case = TRUE), "Phf20", "Scr")
group <- factor(group, levels = c("Scr", "Phf20"))

sample_metadata <- data.frame(
  sample = sample_cols,
  group = group,
  stringsAsFactors = FALSE
)

count_mat <- as.matrix(counts_df[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- make.unique(counts_df$feature_id)

y <- DGEList(counts = count_mat, group = group)
keep <- filterByExpr(y, group = group)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")

design <- model.matrix(~ group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "groupPhf20")

res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
res$feature_id <- counts_df$feature_id[keep][match(rownames(res), rownames(y$counts))]
res$feature_type <- counts_df$feature_type[keep][match(rownames(res), rownames(y$counts))]
res$mean_count <- rowMeans(count_mat[keep, , drop = FALSE])[match(rownames(res), rownames(y$counts))]
res <- res[, c("feature_id", "feature_type", "logFC", "logCPM", "F", "PValue", "FDR", "mean_count")]

res$regulation <- "Not significant"
res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC >= 1] <- "Up in Phf20"
res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC <= -1] <- "Down in Phf20"

gene_res <- res[res$feature_type == "gene", , drop = FALSE]
te_res <- res[res$feature_type == "TE", , drop = FALSE]

write.csv(sample_metadata, file.path(files_dir, paste0(prefix, "_edgeR_sample_metadata.csv")), row.names = FALSE, quote = FALSE)
write.csv(res, file.path(files_dir, paste0(prefix, "_edgeR_DE_gene_TE.csv")), row.names = FALSE, quote = FALSE)
write.csv(gene_res, file.path(files_dir, paste0(prefix, "_edgeR_DE_gene.csv")), row.names = FALSE, quote = FALSE)
write.csv(te_res, file.path(files_dir, paste0(prefix, "_edgeR_DE_TE.csv")), row.names = FALSE, quote = FALSE)

summary_df <- data.frame(
  set = c("gene_TE", "gene", "TE"),
  tested_features = c(nrow(res), nrow(gene_res), nrow(te_res)),
  up_in_Phf20_FDR_0.05_logFC_1 = c(
    sum(res$regulation == "Up in Phf20"),
    sum(gene_res$regulation == "Up in Phf20"),
    sum(te_res$regulation == "Up in Phf20")
  ),
  down_in_Phf20_FDR_0.05_logFC_1 = c(
    sum(res$regulation == "Down in Phf20"),
    sum(gene_res$regulation == "Down in Phf20"),
    sum(te_res$regulation == "Down in Phf20")
  ),
  stringsAsFactors = FALSE
)
write.csv(summary_df, file.path(files_dir, paste0(prefix, "_edgeR_DE_summary.csv")), row.names = FALSE, quote = FALSE)

plot_volcano <- function(df, label, out_prefix) {
  df <- df[is.finite(df$logFC) & !is.na(df$PValue), , drop = FALSE]
  df$neg_log10_fdr <- -log10(pmax(df$FDR, .Machine$double.xmin))
  df$label <- ""
  label_idx <- order(df$FDR, -abs(df$logFC))[seq_len(min(15, nrow(df)))]
  df$label[label_idx] <- df$feature_id[label_idx]

  p <- ggplot(df, aes(x = logFC, y = neg_log10_fdr, color = regulation)) +
    geom_point(alpha = 0.65, size = 1.1) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.3, color = "grey40") +
    geom_text_repel(aes(label = label), max.overlaps = Inf, size = 3, min.segment.length = 0) +
    scale_color_manual(
      values = c(
        "Up in Phf20" = "#D73027",
        "Down in Phf20" = "#4575B4",
        "Not significant" = "grey70"
      ),
      breaks = c("Up in Phf20", "Down in Phf20", "Not significant")
    ) +
    labs(
      title = paste0(label, " differential expression"),
      subtitle = "edgeR QL GLM; logFC = Phf20 / Scr",
      x = "log2 fold change",
      y = "-log10(FDR)",
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )

  ggsave(file.path(plots_dir, paste0(out_prefix, ".pdf")), p, width = 10, height = 6)
  ggsave(file.path(plots_dir, paste0(out_prefix, ".png")), p, width = 10, height = 6, dpi = 300)
}

plot_volcano(res, "Gene + TE", paste0(prefix, "_edgeR_volcano_gene_TE"))
plot_volcano(gene_res, "Gene", paste0(prefix, "_edgeR_volcano_gene"))
plot_volcano(te_res, "TE", paste0(prefix, "_edgeR_volcano_TE"))

message("Input samples:")
print(sample_metadata, row.names = FALSE)
message("Features before filter: ", nrow(counts_df))
message("Features tested after filterByExpr: ", nrow(res))
message("Output files directory: ", files_dir)
message("Output plots directory: ", plots_dir)
print(summary_df, row.names = FALSE)
