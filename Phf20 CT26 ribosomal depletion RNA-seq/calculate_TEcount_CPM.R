#!/usr/bin/env Rscript

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
files_dir <- file.path(workdir, "files")
output_prefix <- "Phf20_ribo_deleption"

count_file <- file.path(files_dir, paste0(output_prefix, "_TEcount_merged_gene_TE_counts.csv"))
if (!file.exists(count_file)) {
  stop("Merged count file not found: ", count_file)
}

counts <- read.csv(
  count_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("feature_id", "feature_type") %in% colnames(counts))) {
  stop("Input matrix must contain feature_id and feature_type columns.")
}

sample_cols <- setdiff(colnames(counts), c("feature_id", "feature_type"))
count_mat <- as.matrix(counts[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"

library_sizes <- colSums(count_mat, na.rm = TRUE)
if (any(library_sizes <= 0)) {
  stop("Found sample(s) with non-positive library size: ",
       paste(names(library_sizes)[library_sizes <= 0], collapse = ", "))
}

cpm_mat <- sweep(count_mat, 2, library_sizes, FUN = "/") * 1e6
cpm <- data.frame(
  feature_id = counts$feature_id,
  feature_type = counts$feature_type,
  cpm_mat,
  check.names = FALSE
)

library_df <- data.frame(
  sample = names(library_sizes),
  library_size = as.numeric(library_sizes),
  cpm_column_sum = as.numeric(colSums(cpm_mat, na.rm = TRUE)),
  stringsAsFactors = FALSE
)

write.csv(
  cpm,
  file = file.path(files_dir, paste0(output_prefix, "_TEcount_merged_gene_TE_CPM.csv")),
  quote = FALSE,
  row.names = FALSE
)

message("Input features: ", nrow(counts))
message("Samples: ", length(sample_cols))
message("Gene features: ", sum(counts$feature_type == "gene"))
message("TE features: ", sum(counts$feature_type == "TE"))
message("Library sizes:")
print(library_df, row.names = FALSE)
message("Output directory: ", files_dir)
