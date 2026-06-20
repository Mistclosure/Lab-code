suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(Seurat)
  library(CRDscore)
  library(qs)
  library(dplyr)
  library(data.table)
  library(Matrix)
})

set.seed(1234)

project_name <- "GSE178318_GSE231559_GSE225857"
work_dir <- file.path("/mnt/disk1/qiuzerui/downloads/CRC", project_name)
input_qs <- file.path(work_dir, "qs", "GSE178318_GSE231559_GSE225857_malignant_primary_RNA_harmony.qs")
signature_file <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table2_primary_cilia_single_gene_confirmed.csv"

files_dir <- file.path(work_dir, "files", "CRDscore")
plots_dir <- file.path(work_dir, "plots", "CRDscore")
qs_dir <- file.path(work_dir, "qs")
logs_dir <- file.path(work_dir, "logs")

dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "02_CRDscore_merged_primary_source.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

signature_name <- sub("\\.csv$", "", basename(signature_file))
signature_prefix <- gsub("[^A-Za-z0-9_]+", "_", signature_name)

message("Step 1: reading merged malignant primary Seurat object")
if (!file.exists(input_qs)) {
  stop("Input qs not found: ", input_qs)
}
obj <- qread(input_qs)
DefaultAssay(obj) <- "RNA"

required_meta <- c("source_dataset", "source_sample", "source_sample_raw", "patient_id",
                   "patient_id_raw", "metastasis_status", "primary_label")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

message("Object: ", nrow(obj), " genes x ", ncol(obj), " cells")
print(table(obj$source_dataset, obj$primary_label, useNA = "ifany"))
print(table(obj$source_sample, useNA = "ifany"))

message("Step 2: NormalizeData using logCPM scale.factor = 1e6")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1000000)

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts) > 0
if (sum(keep_genes) < nrow(obj)) {
  message("Filtering all-zero genes: ", nrow(obj) - sum(keep_genes), " removed")
  obj <- obj[keep_genes, ]
}
rm(counts)
gc()

message("Step 3: preparing log2 expression matrix")
expr_log2 <- LayerData(obj, assay = "RNA", layer = "data") / log(2)

message("Step 4: reading signature genes")
signature_df <- fread(signature_file, data.table = FALSE)
gene_col <- if ("Gene" %in% colnames(signature_df)) "Gene" else colnames(signature_df)[1]
target_genes_raw <- unique(as.character(signature_df[[gene_col]]))
target_genes <- intersect(target_genes_raw, rownames(obj))

if (length(target_genes) <= 5) {
  stop("Too few overlapping signature genes for CRDscore: ", length(target_genes))
}

message("Signature genes in file: ", length(target_genes_raw))
message("Signature genes matched in object: ", length(target_genes))

message("Step 5: optimizing CRDscore parameters for merged data")
requested_n_bins <- 50L
requested_num_rounds <- 1000L

available_genes <- nrow(obj)
unique_gene_means <- length(unique(round(Matrix::rowMeans(expr_log2), digits = 8)))
n_bins <- min(requested_n_bins, max(5L, floor(available_genes / 100L)), max(5L, unique_gene_means - 1L))
num_rounds <- requested_num_rounds

parameter_summary <- data.frame(
  input_qs = input_qs,
  signature_file = signature_file,
  signature_name = signature_name,
  cells = ncol(obj),
  genes_after_filter = nrow(obj),
  signature_genes_raw = length(target_genes_raw),
  signature_genes_matched = length(target_genes),
  requested_n_bins = requested_n_bins,
  final_n_bins = n_bins,
  requested_num_rounds = requested_num_rounds,
  final_num_rounds = num_rounds,
  group_primary_label = paste(levels(factor(obj$primary_label)), collapse = ";"),
  group_source_sample = paste(levels(factor(obj$source_sample)), collapse = ";")
)
write.csv(parameter_summary, file.path(files_dir, "CRDscore_parameter_summary.csv"), row.names = FALSE)

message("Using n.bins = ", n_bins, ", num.rounds = ", num_rounds)

message("Step 6: calculating CRDscore")
score <- cal_CRDscore(
  expr = expr_log2,
  n.bins = n_bins,
  num.rounds = num_rounds,
  circadians = target_genes,
  study.type = "scRNAseq",
  seed = TRUE
)
gc()

score_df <- data.frame(
  cell_id = names(score),
  CRDscore = as.numeric(score),
  row.names = NULL
)

meta_df <- obj@meta.data %>%
  tibble::rownames_to_column("cell_id")

score_meta <- score_df %>%
  left_join(meta_df, by = "cell_id")

if (any(is.na(score_meta$primary_label))) {
  stop("Some scored cells failed to match metadata")
}

write.csv(
  score_meta,
  file.path(files_dir, paste0(signature_prefix, "_merged_primary_source_CRDscore.csv")),
  row.names = FALSE,
  quote = FALSE
)

obj$CRDscore <- score_df$CRDscore[match(colnames(obj), score_df$cell_id)]
qsave(obj, file.path(qs_dir, paste0(project_name, "_malignant_primary_RNA_harmony_CRDscore.qs")))

message("Step 7: summarizing CRDscore by primary source groups")
summary_by_primary_label <- score_meta %>%
  group_by(primary_label) %>%
  summarise(
    cells = n(),
    patients = n_distinct(patient_id),
    samples = n_distinct(source_sample),
    mean_CRDscore = mean(CRDscore, na.rm = TRUE),
    median_CRDscore = median(CRDscore, na.rm = TRUE),
    sd_CRDscore = sd(CRDscore, na.rm = TRUE),
    .groups = "drop"
  )

summary_by_source_sample <- score_meta %>%
  group_by(source_dataset, source_sample, source_sample_raw, patient_id, patient_id_raw,
           metastasis_status, primary_label) %>%
  summarise(
    cells = n(),
    mean_CRDscore = mean(CRDscore, na.rm = TRUE),
    median_CRDscore = median(CRDscore, na.rm = TRUE),
    sd_CRDscore = sd(CRDscore, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(summary_by_primary_label, file.path(files_dir, "CRDscore_summary_by_primary_label.csv"), row.names = FALSE)
write.csv(summary_by_source_sample, file.path(files_dir, "CRDscore_summary_by_source_sample.csv"), row.names = FALSE)

wilcox_by_primary_label <- score_meta %>%
  select(CRDscore, primary_label) %>%
  filter(!is.na(primary_label))

if (n_distinct(wilcox_by_primary_label$primary_label) == 2) {
  wilcox_res <- wilcox.test(CRDscore ~ primary_label, data = wilcox_by_primary_label)
  wilcox_table <- data.frame(
    group1 = levels(factor(wilcox_by_primary_label$primary_label))[1],
    group2 = levels(factor(wilcox_by_primary_label$primary_label))[2],
    method = "wilcox.test",
    p_value = wilcox_res$p.value
  )
  write.csv(wilcox_table, file.path(files_dir, "CRDscore_wilcox_primary_label.csv"), row.names = FALSE)
}

plot_box <- function(data, group_col, title, output_prefix, width = 7, height = 6) {
  plot_data <- data %>%
    transmute(
      CRDscore = CRDscore,
      Type = .data[[group_col]]
    ) %>%
    filter(!is.na(Type), !is.na(CRDscore))

  group_levels <- levels(factor(plot_data$Type))
  plot_data$Type <- factor(plot_data$Type, levels = group_levels)

  comparisons <- list()
  if (length(group_levels) >= 2 && length(group_levels) <= 8) {
    comp <- combn(group_levels, 2)
    comparisons <- lapply(seq_len(ncol(comp)), function(i) comp[, i])
  }

  p <- ggplot(plot_data, aes(x = Type, y = CRDscore, color = Type)) +
    stat_boxplot(geom = "errorbar", width = 0.6) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.7, width = 0.7, fatten = 0.7) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    ggtitle(title) +
    ylab("CRDScore") +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 12, face = "bold", color = "black"),
      axis.text.y = element_text(size = 13, face = "bold", color = "black"),
      axis.title.y = element_text(size = 14, face = "bold", color = "black"),
      plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
    )

  if (length(comparisons) > 0) {
    p <- p + stat_compare_means(comparisons = comparisons, method = "wilcox.test")
  }

  ggsave(file.path(plots_dir, paste0(output_prefix, ".png")), plot = p, width = width, height = height, dpi = 300, bg = "white")
  ggsave(file.path(plots_dir, paste0(output_prefix, ".pdf")), plot = p, width = width, height = height)
  p
}

p_primary_label <- plot_box(
  score_meta,
  "primary_label",
  paste0(signature_name, " CRDscore by primary label"),
  paste0(signature_prefix, "_CRDscore_by_primary_label"),
  width = 7,
  height = 6
)

p_source_sample <- plot_box(
  score_meta,
  "source_sample",
  paste0(signature_name, " CRDscore by primary source"),
  paste0(signature_prefix, "_CRDscore_by_source_sample"),
  width = 10,
  height = 6
)

p_dataset <- plot_box(
  score_meta,
  "source_dataset",
  paste0(signature_name, " CRDscore by dataset"),
  paste0(signature_prefix, "_CRDscore_by_dataset"),
  width = 7,
  height = 6
)

message("Done.")
