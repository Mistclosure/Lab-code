# ==============================================================================
# LUAD: read selected raw files, create Seurat objects, QC, doublet removal, save qs
# Reference: /home/zerui/code/CRC/GSE132465/R/01_GSE132465_data_and_qc.r
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(Matrix)
  library(data.table)
  library(qs)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(1234)

base_dir <- "/mnt/disk1/qiuzerui/downloads/LUAD"
project_dir <- "/home/zerui/code/LUAD"

gse131907_dir <- file.path(base_dir, "GSE131907")
gse123904_dir <- file.path(base_dir, "GSE123904")

clinical_file <- file.path(base_dir, "LUAD_clinical_merged_from_image.csv")
gse131907_counts_file <- file.path(gse131907_dir, "GSE131907_Lung_Cancer_raw_UMI_matrix.txt")
gse131907_annotation_file <- file.path(gse131907_dir, "GSE131907_Lung_Cancer_cell_annotation.txt")

files_dir <- file.path(base_dir, "files")
plots_dir <- file.path(base_dir, "plots")
qs_dir <- file.path(base_dir, "qs", "Seurat")
logs_dir <- file.path(base_dir, "logs")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "01_LUAD_data_and_qc.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("Step 0: loading sample metadata")
clinical <- fread(clinical_file, data.table = FALSE, na.strings = c("", "NA"))

required_cols <- c(
  "Patients ID in this study",
  "Patients ID in GSE131907 / GSE123904",
  "Source_Dataset"
)
missing_cols <- setdiff(required_cols, colnames(clinical))
if (length(missing_cols) > 0) {
  stop("Missing required clinical columns: ", paste(missing_cols, collapse = ", "))
}

clinical <- clinical %>%
  mutate(
    study_id = .data[["Patients ID in this study"]],
    source_id = .data[["Patients ID in GSE131907 / GSE123904"]],
    dataset = .data[["Source_Dataset"]],
    sample_group = case_when(
      str_starts(study_id, "P") ~ "Primary",
      str_starts(study_id, "M") ~ "Metastasis",
      str_starts(study_id, "C") ~ "Control",
      TRUE ~ NA_character_
    )
  )

write.csv(clinical, file.path(files_dir, "LUAD_selected_sample_metadata.csv"), row.names = FALSE)

add_qc_metrics <- function(obj) {
  mt_pattern <- if (any(grepl("^MT-", rownames(obj)))) "^MT-" else "^mt-"
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = mt_pattern)
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^[Rr][Pp][SsLl]")
  obj
}

basic_qc_and_remove_ribo <- function(obj, min_features = 200, max_features = 6000, max_mt = 15) {
  obj <- add_qc_metrics(obj)
  obj <- subset(
    obj,
    subset = nFeature_RNA >= min_features &
      nFeature_RNA <= max_features &
      percent.mt < max_mt
  )

  ribo_genes <- grep("^[Rr][Pp][SsLl]", rownames(obj), value = TRUE)
  if (length(ribo_genes) > 0) {
    obj <- subset(obj, features = setdiff(rownames(obj), ribo_genes))
  }
  obj
}

run_scdblfinder_by_sample <- function(obj, split_by = "study_id") {
  if (!split_by %in% colnames(obj@meta.data)) {
    stop("Metadata does not contain split column: ", split_by)
  }

  obj_list <- SplitObject(obj, split.by = split_by)
  rm(obj)
  gc()

  clean_list <- lapply(names(obj_list), function(sample_name) {
    sample_obj <- obj_list[[sample_name]]
    message("scDblFinder sample: ", sample_name, " | cells after basic QC: ", ncol(sample_obj))

    if (ncol(sample_obj) < 50) {
      message("Skip scDblFinder because cell number < 50: ", sample_name)
      sample_obj$scDblFinder_class <- "not_run_lt50cells"
      return(sample_obj)
    }

    sample_obj <- JoinLayers(sample_obj)
    sample_obj <- NormalizeData(sample_obj, verbose = FALSE)
    sample_obj <- FindVariableFeatures(sample_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    sample_obj <- ScaleData(sample_obj, verbose = FALSE)
    sample_obj <- RunPCA(sample_obj, verbose = FALSE)
    sample_obj <- FindNeighbors(sample_obj, dims = 1:20, verbose = FALSE)
    sample_obj <- FindClusters(sample_obj, resolution = 0.5, verbose = FALSE)

    tryCatch({
      sce <- as.SingleCellExperiment(sample_obj)
      sce <- scDblFinder(sce, clusters = TRUE)
      sample_obj$scDblFinder_class <- sce$scDblFinder.class
      n_doublet <- sum(sample_obj$scDblFinder_class == "doublet", na.rm = TRUE)
      message(
        "scDblFinder doublets: ", sample_name, " = ", n_doublet,
        " (", round(n_doublet / ncol(sample_obj) * 100, 2), "%)"
      )
      sample_obj <- subset(sample_obj, subset = scDblFinder_class == "singlet")
    }, error = function(e) {
      message("scDblFinder failed for ", sample_name, ": ", e$message)
      sample_obj$scDblFinder_class <- "not_run_error"
    })

    qsave(sample_obj, file.path(qs_dir, paste0("LUAD_", sample_name, "_qc.qs")))
    gc()
    sample_obj
  })

  clean_list[!vapply(clean_list, is.null, logical(1))]
}

merge_seurat_list <- function(obj_list, project_name) {
  if (length(obj_list) == 0) {
    return(NULL)
  }
  if (length(obj_list) == 1) {
    return(JoinLayers(obj_list[[1]]))
  }
  merged <- merge(obj_list[[1]], y = obj_list[-1], project = project_name)
  JoinLayers(merged)
}

make_unique_cells <- function(cells, prefix) {
  paste(prefix, make.unique(as.character(cells)), sep = "_")
}

repair_gse123904_metadata <- function(obj, clinical) {
  if (!"orig.ident" %in% colnames(obj@meta.data)) {
    return(obj)
  }

  md <- obj@meta.data
  for (col in c("source_dataset", "source_id", "study_id", "dataset", "sample_group")) {
    if (!col %in% colnames(md)) {
      md[[col]] <- NA_character_
    }
  }

  missing_gse123904 <- md$orig.ident == "GSE123904" &
    (is.na(md$study_id) | is.na(md$source_id) | is.na(md$source_dataset))

  if (!any(missing_gse123904, na.rm = TRUE)) {
    obj@meta.data <- md
    return(obj)
  }

  extracted_source_id <- str_match(rownames(md)[missing_gse123904], "^GSE123904_([^_]+)_")[, 2]
  clinical_match <- match(extracted_source_id, clinical$source_id)

  md$source_dataset[missing_gse123904] <- "GSE123904"
  md$source_id[missing_gse123904] <- extracted_source_id
  md$dataset[missing_gse123904] <- "GSE123904"
  md$study_id[missing_gse123904] <- clinical$study_id[clinical_match]
  md$sample_group[missing_gse123904] <- clinical$sample_group[clinical_match]

  obj@meta.data <- md
  obj
}

add_clinical_metadata <- function(obj, clinical) {
  md <- obj@meta.data
  old_clinical_cols <- grep("^clinical__", colnames(md), value = TRUE)
  if (length(old_clinical_cols) > 0) {
    md <- md[, setdiff(colnames(md), old_clinical_cols), drop = FALSE]
  }

  clinical_meta <- clinical %>%
    distinct(study_id, source_id, .keep_all = TRUE) %>%
    rename_with(
      ~ paste0("clinical__", make.names(.x)),
      .cols = -c(study_id, source_id)
    )

  md <- md %>%
    tibble::rownames_to_column("cell_barcode") %>%
    left_join(clinical_meta, by = c("study_id", "source_id")) %>%
    as.data.frame()
  rownames(md) <- md$cell_barcode
  md$cell_barcode <- NULL

  obj@meta.data <- md
  obj
}

read_gse131907_selected <- function() {
  message("Step 1: reading selected GSE131907 cells")
  selected <- clinical %>%
    filter(dataset == "GSE131907") %>%
    pull(source_id) %>%
    unique()

  if (length(selected) == 0) {
    message("No GSE131907 samples selected")
    return(NULL)
  }

  annotation <- fread(gse131907_annotation_file, data.table = FALSE)
  rownames(annotation) <- annotation$Index
  annotation <- annotation %>%
    filter(Sample %in% selected)

  if (nrow(annotation) == 0) {
    stop("No selected GSE131907 cells found in annotation")
  }

  selected_cells <- annotation$Index
  matrix_header <- names(fread(gse131907_counts_file, nrows = 0))
  selected_cols <- c("Index", intersect(selected_cells, matrix_header))
  missing_cells <- setdiff(selected_cells, matrix_header)
  if (length(missing_cells) > 0) {
    warning("GSE131907 cells absent from matrix: ", length(missing_cells))
  }

  counts_dt <- fread(
    gse131907_counts_file,
    select = selected_cols,
    data.table = FALSE,
    showProgress = TRUE
  )
  rownames(counts_dt) <- counts_dt$Index
  counts_dt$Index <- NULL

  counts <- as(as.matrix(counts_dt), "dgCMatrix")
  rownames(counts) <- make.unique(rownames(counts))
  rm(counts_dt)
  gc()

  annotation <- annotation[colnames(counts), , drop = FALSE]
  annotation <- annotation %>%
    mutate(
      source_dataset = "GSE131907",
      source_id = Sample,
      original_cell_id = rownames(annotation)
    ) %>%
    left_join(
      clinical,
      by = c("source_id" = "source_id"),
      suffix = c("", ".clinical")
    )
  rownames(annotation) <- annotation$original_cell_id

  colnames(counts) <- make_unique_cells(colnames(counts), "GSE131907")
  rownames(annotation) <- colnames(counts)

  obj <- CreateSeuratObject(counts = counts, meta.data = annotation, project = "LUAD_GSE131907")
  rm(counts, annotation)
  gc()

  qsave(obj, file.path(qs_dir, "LUAD_GSE131907_raw_selected.qs"))
  obj
}

pick_gse123904_file <- function(source_id, sample_group) {
  pattern <- paste0("MSK_", source_id, "_.*dense\\.csv\\.gz$")
  candidates <- list.files(gse123904_dir, pattern = pattern, full.names = TRUE)
  if (length(candidates) == 0) {
    candidates <- list.files(gse123904_dir, pattern = paste0(source_id, ".*dense\\.csv\\.gz$"), full.names = TRUE)
  }
  if (length(candidates) == 0) {
    stop("No GSE123904 dense file found for ", source_id)
  }

  if (sample_group == "Primary") {
    hit <- candidates[str_detect(basename(candidates), "PRIMARY_TUMOUR")]
    if (length(hit) > 0) return(hit[1])
  }
  if (sample_group %in% c("Metastasis", "Control")) {
    hit <- candidates[str_detect(basename(candidates), "METASTASIS|NORMAL")]
    if (length(hit) > 0) return(hit[1])
  }

  candidates[1]
}

read_gse123904_selected <- function() {
  message("Step 2: reading selected GSE123904 dense matrices")
  selected_meta <- clinical %>%
    filter(dataset == "GSE123904") %>%
    distinct(study_id, source_id, sample_group, .keep_all = TRUE)

  if (nrow(selected_meta) == 0) {
    message("No GSE123904 samples selected")
    return(NULL)
  }

  obj_list <- lapply(seq_len(nrow(selected_meta)), function(i) {
    sample_meta <- selected_meta[i, , drop = FALSE]
    source_id <- sample_meta$source_id
    sample_group <- sample_meta$sample_group
    input_file <- pick_gse123904_file(source_id, sample_group)

    message("Reading GSE123904 sample: ", sample_meta$study_id, " / ", source_id, " from ", basename(input_file))
    dense <- fread(input_file, data.table = FALSE, check.names = FALSE)

    cell_ids <- dense[[1]]
    dense[[1]] <- NULL
    counts <- t(as.matrix(dense))
    rownames(counts) <- make.unique(colnames(dense))
    colnames(counts) <- make_unique_cells(cell_ids, paste0("GSE123904_", source_id))
    counts <- as(counts, "dgCMatrix")

    meta <- as.data.frame(sample_meta[rep(1, ncol(counts)), , drop = FALSE])
    rownames(meta) <- colnames(counts)
    meta$source_dataset <- "GSE123904"
    meta$original_cell_id <- cell_ids
    meta$source_file <- basename(input_file)

    obj <- CreateSeuratObject(counts = counts, meta.data = meta, project = "LUAD_GSE123904")
    qsave(obj, file.path(qs_dir, paste0("LUAD_GSE123904_", sample_meta$study_id, "_raw.qs")))
    rm(dense, counts, meta)
    gc()
    obj
  })

  merged <- merge_seurat_list(obj_list, "LUAD_GSE123904")
  rm(obj_list)
  gc()

  qsave(merged, file.path(qs_dir, "LUAD_GSE123904_raw_selected.qs"))
  merged
}

raw_merged_qs <- file.path(qs_dir, "LUAD_selected_raw_merged.qs")
if (file.exists(raw_merged_qs)) {
  message("Step 3: loading existing selected raw merged qs")
  sc_raw <- qread(raw_merged_qs)
} else {
  message("Step 3: creating selected raw Seurat objects")
  gse131907_obj <- read_gse131907_selected()
  gse123904_obj <- read_gse123904_selected()

  raw_list <- list(gse131907_obj, gse123904_obj)
  raw_list <- raw_list[!vapply(raw_list, is.null, logical(1))]
  if (length(raw_list) == 0) {
    stop("No Seurat objects were created")
  }

  sc_raw <- merge_seurat_list(raw_list, "LUAD_selected_raw")
  qsave(sc_raw, raw_merged_qs)
  rm(raw_list, gse131907_obj, gse123904_obj)
  gc()
}

sc_raw <- repair_gse123904_metadata(sc_raw, clinical)
sc_raw <- add_clinical_metadata(sc_raw, clinical)
qsave(sc_raw, raw_merged_qs)

qc_before <- sc_raw@meta.data %>%
  dplyr::count(source_dataset, study_id, source_id, name = "cells_raw")
write.csv(qc_before, file.path(files_dir, "LUAD_cells_raw_by_sample.csv"), row.names = FALSE)

message("Step 4: basic QC and ribosomal gene removal")
sc_qc <- basic_qc_and_remove_ribo(sc_raw)
sc_qc <- add_clinical_metadata(sc_qc, clinical)
rm(sc_raw)
gc()
qsave(sc_qc, file.path(qs_dir, "LUAD_selected_basic_qc_no_ribo.qs"))

qc_after_basic <- sc_qc@meta.data %>%
  dplyr::count(source_dataset, study_id, source_id, name = "cells_after_basic_qc")
write.csv(qc_after_basic, file.path(files_dir, "LUAD_cells_after_basic_qc_by_sample.csv"), row.names = FALSE)

message("Step 5: scDblFinder by study sample")
sc_clean_list <- run_scdblfinder_by_sample(sc_qc, split_by = "study_id")
rm(sc_qc)
gc()

message("Step 6: final merge and qs save")
sc_final <- merge_seurat_list(sc_clean_list, "LUAD_selected_qc_cleaned")
sc_final <- add_clinical_metadata(sc_final, clinical)
qsave(sc_final, file.path(qs_dir, "LUAD_selected_qc_cleaned.qs"))

qc_final <- sc_final@meta.data %>%
  dplyr::count(source_dataset, study_id, source_id, sample_group, scDblFinder_class, name = "cells_final")
write.csv(qc_final, file.path(files_dir, "LUAD_cells_final_by_sample.csv"), row.names = FALSE)

message("Final cells: ", ncol(sc_final), " | genes: ", nrow(sc_final))
print(table(sc_final$study_id))

rm(sc_clean_list)
gc()

message("Done. Final qs: ", file.path(qs_dir, "LUAD_selected_qc_cleaned.qs"))
