#!/usr/bin/env Rscript

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20 CT26 ribosomal depletion RNA-seq"
analysis_dir <- file.path(workdir, "sense_antisense_ERV")
results_dir <- file.path(analysis_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(workdir, "files", "Phf20_ribo_deleption_edgeR_sample_metadata.csv")
sense_dir <- file.path(analysis_dir, "counts", "TEcount_sense")
antisense_dir <- file.path(analysis_dir, "counts", "TEcount_antisense")
target_classes <- c("LINE", "SINE", "LTR")

if (!file.exists(metadata_file)) {
  stop("Metadata not found: ", metadata_file, "\nPlease confirm sample metadata before merging.")
}
metadata <- read.csv(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("sample", "group") %in% colnames(metadata))) {
  stop("Metadata must contain columns: sample, group")
}

read_one <- function(path, sample_name) {
  dat <- read.delim(path, header = TRUE, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
  if (ncol(dat) != 2) stop("Unexpected column number in ", path)
  colnames(dat) <- c("feature_id", sample_name)
  dat[[sample_name]] <- as.numeric(dat[[sample_name]])
  dat
}

merge_dir <- function(count_dir, metadata) {
  files <- list.files(count_dir, pattern = "\\.cntTable$", full.names = TRUE)
  if (length(files) == 0) stop("No .cntTable files found in: ", count_dir)
  names(files) <- make.names(sub("\\.cntTable$", "", basename(files)), unique = TRUE)

  missing <- setdiff(metadata$sample, names(files))
  extra <- setdiff(names(files), metadata$sample)
  if (length(missing) > 0) stop("Missing samples in ", count_dir, ": ", paste(missing, collapse = ", "))
  if (length(extra) > 0) message("Ignoring extra samples in ", count_dir, ": ", paste(extra, collapse = ", "))

  files <- files[metadata$sample]
  count_list <- Map(read_one, files, metadata$sample)
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), count_list)
  merged[is.na(merged)] <- 0
  for (sample in metadata$sample) merged[[sample]] <- as.numeric(merged[[sample]])
  merged
}

filter_te <- function(dat) {
  dat[!grepl("^ENSMUSG", dat$feature_id), , drop = FALSE]
}

add_te_class <- function(dat) {
  dat$TE_class <- sub("^.*:([^:]+)$", "\\1", dat$feature_id)
  dat
}

filter_line_sine_ltr <- function(dat) {
  dat <- add_te_class(dat)
  dat[dat$TE_class %in% target_classes, c("feature_id", "TE_class", metadata$sample), drop = FALSE]
}

sense_all <- merge_dir(sense_dir, metadata)
antisense_all <- merge_dir(antisense_dir, metadata)

sense_te <- filter_te(sense_all)
antisense_te <- filter_te(antisense_all)
sense_target <- filter_line_sine_ltr(sense_te)
antisense_target <- filter_line_sine_ltr(antisense_te)

write.table(metadata, file.path(results_dir, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sense_te, file.path(results_dir, "TE_sense_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(antisense_te, file.path(results_dir, "TE_antisense_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(sense_target, file.path(results_dir, "TE_sense_LINE_SINE_LTR_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(antisense_target, file.path(results_dir, "TE_antisense_LINE_SINE_LTR_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Metadata samples: ", paste(metadata$sample, collapse = ", "))
message("Sense TE rows: ", nrow(sense_te), "; sense LINE/SINE/LTR rows: ", nrow(sense_target))
message("Antisense TE rows: ", nrow(antisense_te), "; antisense LINE/SINE/LTR rows: ", nrow(antisense_target))
message("Output directory: ", results_dir)
