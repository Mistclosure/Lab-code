#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
})

workdir <- "/mnt/disk1/qiuzerui/expriments/Phf20_1.23"
analysis_dir <- file.path(workdir, "sense_antisense_ERV")
results_dir <- file.path(analysis_dir, "results")
files_dir <- file.path(workdir, "files")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(results_dir, "metadata.tsv")
fallback_metadata_file <- file.path(files_dir, "Phf20_1.23_edgeR_sample_metadata.csv")
sense_dir <- file.path(analysis_dir, "counts", "TEcount_sense")
antisense_dir <- file.path(analysis_dir, "counts", "TEcount_antisense")

infer_metadata <- function(count_dir) {
  files <- list.files(count_dir, pattern = "\\.cntTable$", full.names = FALSE)
  if (length(files) == 0) stop("No .cntTable files found for metadata inference: ", count_dir)
  raw_sample <- sub("\\.cntTable$", "", files)
  sample <- make.names(raw_sample, unique = TRUE)
  group <- ifelse(grepl("Phf20", sample, ignore.case = TRUE), "Phf20",
                  ifelse(grepl("Scr", sample, ignore.case = TRUE), "Scr", NA))
  if (any(is.na(group))) {
    stop("Cannot infer group for samples: ", paste(sample[is.na(group)], collapse = ", "))
  }
  metadata <- data.frame(sample = sample, group = group, stringsAsFactors = FALSE)
  metadata[order(factor(metadata$group, levels = c("Scr", "Phf20")), metadata$sample), ]
}

if (file.exists(metadata_file)) {
  metadata <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
} else if (file.exists(fallback_metadata_file)) {
  metadata <- read.csv(fallback_metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
} else {
  metadata <- infer_metadata(sense_dir)
}
metadata$sample <- make.names(metadata$sample, unique = TRUE)
metadata$group <- factor(metadata$group, levels = c("Scr", "Phf20"))
if (!all(c("sample", "group") %in% colnames(metadata))) {
  stop("Metadata must contain columns: sample and group")
}
if (any(is.na(metadata$group))) {
  stop("Metadata group must be Scr or Phf20")
}
write.table(metadata, file.path(results_dir, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

read_one <- function(path, sample_name) {
  dat <- read.delim(path, header = TRUE, check.names = FALSE, quote = "\"", stringsAsFactors = FALSE)
  if (ncol(dat) != 2) stop("Unexpected column number in ", path)
  colnames(dat) <- c("feature_id", sample_name)
  dat[[sample_name]] <- as.numeric(dat[[sample_name]])
  dat
}

merge_count_dir <- function(count_dir, metadata) {
  files <- list.files(count_dir, pattern = "\\.cntTable$", full.names = TRUE)
  if (length(files) == 0) stop("No .cntTable files found in: ", count_dir)
  names(files) <- make.names(sub("\\.cntTable$", "", basename(files)), unique = TRUE)

  missing <- setdiff(metadata$sample, names(files))
  if (length(missing) > 0) stop("Missing samples in ", count_dir, ": ", paste(missing, collapse = ", "))

  files <- files[metadata$sample]
  count_list <- Map(read_one, files, metadata$sample)
  merged <- Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), count_list)
  merged[is.na(merged)] <- 0
  for (sample in metadata$sample) merged[[sample]] <- as.numeric(merged[[sample]])
  merged
}

annotate_features <- function(feature_id) {
  feature_type <- ifelse(grepl("^ENSMUSG", feature_id), "gene", "TE")
  te_class <- ifelse(feature_type == "TE" & grepl(":", feature_id),
                     sub("^.*:([^:]+)$", "\\1", feature_id), NA)
  data.frame(feature_id = feature_id, feature_type = feature_type, TE_class = te_class,
             stringsAsFactors = FALSE)
}

run_edger <- function(count_dir, label) {
  dat <- merge_count_dir(count_dir, metadata)
  ann <- annotate_features(dat$feature_id)
  count_mat <- as.matrix(dat[, metadata$sample, drop = FALSE])
  rownames(count_mat) <- make.unique(dat$feature_id)
  storage.mode(count_mat) <- "numeric"

  # TEcount multi mode can emit fractional assignments; edgeR accepts numeric counts.
  y <- DGEList(counts = count_mat, group = metadata$group)
  keep <- filterByExpr(y, group = metadata$group)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  design <- model.matrix(~ group, data = metadata)
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "groupPhf20")

  res <- topTags(qlf, n = Inf, sort.by = "PValue")$table
  tested_ann <- ann[keep, , drop = FALSE]
  idx <- match(rownames(res), rownames(y$counts))
  res$feature_id <- tested_ann$feature_id[idx]
  res$feature_type <- tested_ann$feature_type[idx]
  res$TE_class <- tested_ann$TE_class[idx]
  res$mean_count <- rowMeans(count_mat[keep, , drop = FALSE])[idx]
  res <- res[, c("feature_id", "feature_type", "TE_class", "logFC", "logCPM", "F", "PValue", "FDR", "mean_count")]

  res$regulation <- "Not significant"
  res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC >= 1] <- "Up in Phf20"
  res$regulation[!is.na(res$FDR) & res$FDR < 0.05 & res$logFC <= -1] <- "Down in Phf20"

  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  rownames(logcpm) <- res$feature_id[match(rownames(logcpm), rownames(y$counts))]

  out_prefix <- file.path(results_dir, paste0("edgeR_", label, "_all_gene_TE"))
  write.table(res, paste0(out_prefix, "_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(res, paste0(out_prefix, "_results.csv"), row.names = FALSE)
  write.table(
    cbind(feature_id = rownames(logcpm), as.data.frame(logcpm, check.names = FALSE)),
    paste0(out_prefix, "_logCPM_matrix.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  summary_by_type <- aggregate(
    cbind(
      tested_features = rep(1, nrow(res)),
      up_in_Phf20_FDR_0.05_logFC_1 = as.integer(res$regulation == "Up in Phf20"),
      down_in_Phf20_FDR_0.05_logFC_1 = as.integer(res$regulation == "Down in Phf20")
    ) ~ feature_type,
    data = res,
    FUN = sum
  )
  total <- data.frame(
    feature_type = "all_gene_TE",
    tested_features = nrow(res),
    up_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Up in Phf20"),
    down_in_Phf20_FDR_0.05_logFC_1 = sum(res$regulation == "Down in Phf20")
  )
  summary_df <- rbind(total, summary_by_type)
  write.table(summary_df, paste0(out_prefix, "_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

  sig <- res[res$regulation != "Not significant", , drop = FALSE]
  write.table(sig, paste0(out_prefix, "_significant_FDR0.05_logFC1.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)

  message(label, " tested features: ", nrow(res))
  print(summary_df, row.names = FALSE)
  invisible(list(results = res, summary = summary_df))
}

message("Running edgeR for all genes + all TE families")
message("Samples: ", paste(metadata$sample, collapse = ", "))
message("Groups: ", paste(as.character(metadata$group), collapse = ", "))

sense <- run_edger(sense_dir, "sense")
antisense <- run_edger(antisense_dir, "antisense")

message("Done. Output directory: ", results_dir)
