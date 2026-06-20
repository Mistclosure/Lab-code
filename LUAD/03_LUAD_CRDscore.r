# ==============================================================================
# LUAD script 3: CRDscore for CopyKAT-defined malignant cells
# Reference: /home/zerui/code/CRC/GSE132465/R/CRDscore_GSE132465.r
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(tidyverse)
  library(Seurat)
  library(CRDscore)
  library(qs)
  library(data.table)
})

set.seed(1234)

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/LUAD"
setwd(WORK_DIR)

input_qs <- file.path(WORK_DIR, "qs", "Seurat", "LUAD_Malignant_RNA_harmony.qs")
args <- commandArgs(trailingOnly = TRUE)
signature_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/mnt/disk1/qiuzerui/downloads/CRC/signature/table1_primary_cilia_single_gene_confirmed.csv"
}

files_dir <- file.path(WORK_DIR, "files", "clinical_association")
plots_dir <- file.path(WORK_DIR, "plots", "clinical_association")
logs_dir <- file.path(WORK_DIR, "logs")

dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "03_LUAD_CRDscore.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

signature_name <- tools::file_path_sans_ext(basename(signature_file))
signature_file_prefix <- signature_name

theme_crd <- function() {
  theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 15, face = "bold", color = "black"),
      axis.title.y = element_text(size = 15, face = "bold", color = "black"),
      axis.text.x = element_text(size = 15, face = "bold", color = "black", angle = 90, hjust = 1),
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5)
    )
}

make_comparisons <- function(groups) {
  groups <- unique(as.character(groups))
  groups <- groups[!is.na(groups) & groups != ""]
  if (length(groups) < 2) {
    return(list())
  }
  comps <- combn(groups, 2, simplify = FALSE)
  comps
}

plot_score_box <- function(data, group_col, title, output_prefix) {
  plot_data <- data %>%
    transmute(score = .data[["score"]], Type = .data[[group_col]]) %>%
    filter(!is.na(score), !is.na(Type), Type != "")

  if (nrow(plot_data) == 0 || length(unique(plot_data$Type)) < 2) {
    message("Skip plot for ", group_col, ": fewer than two valid groups")
    return(invisible(NULL))
  }

  plot_data$Type <- factor(plot_data$Type, levels = unique(plot_data$Type))
  comparisons <- make_comparisons(levels(plot_data$Type))

  p <- ggplot(plot_data, aes(x = Type, y = score, color = Type)) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test") +
    ggtitle(title) +
    ylab("CRDScore") +
    theme_crd()

  print(p)
  ggsave(file.path(plots_dir, paste0(output_prefix, ".png")), plot = p, width = 6, height = 5, dpi = 300)
  ggsave(file.path(plots_dir, paste0(output_prefix, ".pdf")), plot = p, width = 6, height = 5)
  invisible(p)
}

message("Step 1: reading malignant LUAD Seurat object")
if (!file.exists(input_qs)) {
  stop("Input qs not found: ", input_qs)
}
if (!file.exists(signature_file)) {
  stop("Signature file not found: ", signature_file)
}

pbmc1 <- qread(input_qs)
pbmc1 <- JoinLayers(pbmc1)
DefaultAssay(pbmc1) <- "RNA"

message("Object: ", nrow(pbmc1), " genes x ", ncol(pbmc1), " cells")
message("Assay: ", DefaultAssay(pbmc1))
if ("sample_group" %in% colnames(pbmc1@meta.data)) {
  print(table(pbmc1$sample_group, useNA = "ifany"))
}
if ("clinical__Stage" %in% colnames(pbmc1@meta.data)) {
  print(table(pbmc1$clinical__Stage, useNA = "ifany"))
}

message("Step 2: normalizing data for CRDscore")
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 1000000)
seurat_data <- LayerData(pbmc1, assay = "RNA", layer = "data")
data_log2 <- seurat_data / log(2)

message("Step 3: reading signature genes")
signature_data <- read.csv(signature_file, header = TRUE, check.names = FALSE)
if (!"Gene" %in% colnames(signature_data)) {
  stop("Signature file must contain a Gene column: ", signature_file)
}

target_genes_raw <- unique(as.character(signature_data[["Gene"]]))
target_genes_raw <- target_genes_raw[!is.na(target_genes_raw) & target_genes_raw != ""]
target_genes <- intersect(target_genes_raw, rownames(pbmc1))
missing_genes <- setdiff(target_genes_raw, target_genes)

gene_overlap <- data.frame(
  Gene = target_genes_raw,
  present_in_object = target_genes_raw %in% target_genes
)
write.csv(
  gene_overlap,
  file.path(files_dir, paste0(signature_file_prefix, "_LUAD_CRDscore_gene_overlap.csv")),
  row.names = FALSE,
  quote = FALSE
)

message("Signature genes: ", length(target_genes_raw))
message("Overlapping genes: ", length(target_genes))
if (length(missing_genes) > 0) {
  message("Missing genes: ", paste(missing_genes, collapse = ", "))
}
if (length(target_genes) <= 5) {
  stop("Non-enough-overlapping genes to calculate CRDscore")
}

message("Step 4: calculating CRDscore")
score <- cal_CRDscore(
  expr = data_log2,
  n.bins = 50,
  circadians = target_genes,
  study.type = "scRNAseq"
)
gc()

score <- as.data.frame(score)
score$id <- rownames(score)
colnames(score)[colnames(score) == "score"] <- "score"
if (!"score" %in% colnames(score)) {
  colnames(score)[1] <- "score"
}
score <- score[, c("id", "score")]

message("Step 5: merging CRDscore with metadata")
meta <- pbmc1@meta.data %>%
  tibble::rownames_to_column("id")

rt1 <- score %>%
  left_join(meta, by = "id")

score_file <- file.path(files_dir, paste0(signature_file_prefix, "_LUAD_CRDscore.csv"))
score_only_file <- file.path(files_dir, paste0(signature_file_prefix, "_LUAD_CRDscore_score_only.csv"))
write.csv(rt1, score_file, row.names = FALSE, quote = FALSE)
write.csv(score, score_only_file, row.names = FALSE, quote = FALSE)

message("Step 6: plotting clinical associations")
plot_df <- rt1

if ("clinical__Stage" %in% colnames(plot_df)) {
  plot_df$Stage_for_plot <- plot_df$clinical__Stage
} else if ("Stage" %in% colnames(plot_df)) {
  plot_df$Stage_for_plot <- plot_df$Stage
}

if ("Stage_for_plot" %in% colnames(plot_df)) {
  stage_levels <- c("IA", "IA3", "IB", "IIA", "IIB", "IIIA", "IIIB", "IIIC", "IV")
  present_stage <- intersect(stage_levels, unique(as.character(plot_df$Stage_for_plot)))
  other_stage <- setdiff(unique(as.character(plot_df$Stage_for_plot)), c(present_stage, NA, ""))
  plot_df$Stage_for_plot <- factor(plot_df$Stage_for_plot, levels = c(present_stage, other_stage))
  plot_score_box(
    plot_df,
    "Stage_for_plot",
    paste0(signature_name, "+LUAD+CRDscore (Stage)"),
    paste0(signature_file_prefix, "_LUAD_Stage_plot")
  )
}

if ("sample_group" %in% colnames(plot_df)) {
  group_levels <- c("Control", "Primary", "Metastasis")
  present_group <- intersect(group_levels, unique(as.character(plot_df$sample_group)))
  other_group <- setdiff(unique(as.character(plot_df$sample_group)), c(present_group, NA, ""))
  plot_df$sample_group <- factor(plot_df$sample_group, levels = c(present_group, other_group))
  plot_score_box(
    plot_df,
    "sample_group",
    paste0(signature_name, "+LUAD+CRDscore (Sample group)"),
    paste0(signature_file_prefix, "_LUAD_Sample_group_plot")
  )
}

if ("source_dataset" %in% colnames(plot_df)) {
  plot_score_box(
    plot_df,
    "source_dataset",
    paste0(signature_name, "+LUAD+CRDscore (Dataset)"),
    paste0(signature_file_prefix, "_LUAD_Dataset_plot")
  )
}

message("Done. CRDscore outputs are saved under: ", WORK_DIR)
