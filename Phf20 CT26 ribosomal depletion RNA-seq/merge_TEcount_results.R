#!/usr/bin/env Rscript

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
counts_dir <- file.path(workdir, "counts")
files_dir <- file.path(workdir, "files")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
output_prefix <- "Phf20_ribo_deleption"

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

write.csv(
  merged_with_type,
  file = file.path(files_dir, paste0(output_prefix, "_TEcount_merged_gene_TE_counts.csv")),
  quote = FALSE,
  row.names = FALSE
)

message("Merged files: ", length(cnt_files))
message("Total features: ", nrow(merged_with_type))
message("Gene features: ", sum(merged_with_type$feature_type == "gene"))
message("TE features: ", sum(merged_with_type$feature_type == "TE"))
message("Output directory: ", files_dir)
