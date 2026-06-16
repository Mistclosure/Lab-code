#!/usr/bin/env Rscript

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
counts_dir <- file.path(workdir, "counts")
files_dir <- file.path(workdir, "files")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

cnt_files <- list.files(
  counts_dir,
  pattern = "\\.cntTable$",
  full.names = TRUE
)

if (length(cnt_files) == 0) {
  stop("No .cntTable files found in: ", counts_dir)
}

cnt_files <- sort(cnt_files)

sample_names <- sub("\\.cntTable$", "", basename(cnt_files))
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

gene_matrix <- merged_with_type[merged_with_type$feature_type == "gene", , drop = FALSE]
te_matrix <- merged_with_type[merged_with_type$feature_type == "TE", , drop = FALSE]

metadata <- data.frame(
  sample = sample_names,
  source_file = basename(cnt_files),
  group = ifelse(grepl("Phf20", sample_names, ignore.case = TRUE), "Phf20", "Scr"),
  replicate = sub(".*_(\\d+)_Mixt$", "\\1", sample_names),
  stringsAsFactors = FALSE
)

write.table(
  merged_with_type,
  file = file.path(files_dir, "TEcount_merged_gene_TE_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  merged_with_type,
  file = file.path(files_dir, "TEcount_merged_gene_TE_counts.csv"),
  quote = FALSE,
  row.names = FALSE
)

write.table(
  gene_matrix,
  file = file.path(files_dir, "TEcount_merged_gene_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  te_matrix,
  file = file.path(files_dir, "TEcount_merged_TE_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.csv(
  metadata,
  file = file.path(files_dir, "TEcount_sample_metadata.csv"),
  quote = FALSE,
  row.names = FALSE
)

message("Merged files: ", length(cnt_files))
message("Total features: ", nrow(merged_with_type))
message("Gene features: ", nrow(gene_matrix))
message("TE features: ", nrow(te_matrix))
message("Output directory: ", files_dir)
