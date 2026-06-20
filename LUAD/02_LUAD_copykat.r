# ==============================================================================
# LUAD script 2: RNA workflow + CopyKAT malignant-cell prediction + malignant subset
# Reference: /home/zerui/code/CRC/GSE132465/R/02_GSE132465_copykat.r
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(stringr)
  library(magrittr)
  library(harmony)
  library(qs)
  library(copykat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(1234)

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/LUAD"
setwd(WORK_DIR)

input_qs <- file.path(WORK_DIR, "qs", "Seurat", "LUAD_selected_qc_cleaned.qs")
clinical_file <- file.path(WORK_DIR, "LUAD_clinical_merged_from_image.csv")

files_dir <- file.path(WORK_DIR, "files")
tables_dir <- file.path(files_dir, "tables")
plots_dir <- file.path(WORK_DIR, "plots", "copykat")
seurat_qs_dir <- file.path(WORK_DIR, "qs", "Seurat")
copykat_qs_dir <- file.path(WORK_DIR, "qs", "copykat")
copykat_raw_dir <- file.path(WORK_DIR, "files", "copykat")
logs_dir <- file.path(WORK_DIR, "logs")

dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(seurat_qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(copykat_qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(copykat_raw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "02_LUAD_copykat.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("Step 1: reading QC-cleaned LUAD Seurat object")
if (!file.exists(input_qs)) {
  stop("Input qs not found: ", input_qs)
}
scRNA <- qread(input_qs)
scRNA <- JoinLayers(scRNA)
DefaultAssay(scRNA) <- "RNA"

required_meta <- c("study_id", "source_id", "source_dataset", "sample_group")
missing_meta <- setdiff(required_meta, colnames(scRNA@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

message("Object before CopyKAT: ", nrow(scRNA), " genes x ", ncol(scRNA), " cells")
print(table(scRNA$sample_group, useNA = "ifany"))
print(table(scRNA$study_id, useNA = "ifany"))

# User request: run CopyKAT on all selected LUAD samples, including P/M/C.
message("All selected LUAD samples will be used for CopyKAT")

mt_pattern <- if (any(grepl("^MT-", rownames(scRNA)))) "^MT-" else "^mt-"
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = mt_pattern)

# ------------------------------------------------------------------------------
# 2. Standard RNA workflow and Harmony integration
# ------------------------------------------------------------------------------
message("Step 2: running standard RNA workflow and Harmony integration")
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

batch_var <- "study_id"
scRNA <- RunHarmony(
  scRNA,
  group.by.vars = batch_var,
  max.iter.harmony = 20,
  assay.use = "RNA"
)

pc.num <- 1:30
scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = pc.num)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = pc.num)
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 1.0)

qsave(scRNA, file.path(seurat_qs_dir, "LUAD_all_samples_RNA_harmony_before_copykat.qs"))

p_sample <- DimPlot(scRNA, reduction = "umap", group.by = "study_id", raster = TRUE)
p_cluster <- DimPlot(scRNA, reduction = "umap", label = TRUE, raster = TRUE) + NoLegend()
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_study_id.pdf"), p_sample, width = 8, height = 6)
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_cluster.pdf"), p_cluster, width = 8, height = 6)
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_study_id.png"), p_sample, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_cluster.png"), p_cluster, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3. Run CopyKAT by study sample
# ------------------------------------------------------------------------------
message("Step 3: running CopyKAT by study_id")
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")
sample_ids <- unique(scRNA$study_id)
sample_ids <- sample_ids[!is.na(sample_ids)]

pred_list <- list()
copykat_cores <- 24

for (i in seq_along(sample_ids)) {
  current_id <- sample_ids[i]
  message("========== CopyKAT ", i, "/", length(sample_ids), ": ", current_id, " ==========")

  cell_idx <- which(scRNA$study_id == current_id)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]

  if (ncol(sub_counts) < 50) {
    message("Skip CopyKAT for ", current_id, ": too few cells (", ncol(sub_counts), ")")
    next
  }

  sub_copykat_res <- tryCatch({
    copykat(
      rawmat = as.matrix(sub_counts),
      id.type = "S",
      ngene.chr = 5,
      win.size = 25,
      KS.cut = 0.1,
      sam.name = file.path(copykat_raw_dir, paste0("copykat_", current_id)),
      distance = "euclidean",
      n.cores = copykat_cores
    )
  }, error = function(e) {
    message("CopyKAT failed for ", current_id, ": ", e$message)
    NULL
  })

  if (!is.null(sub_copykat_res) && "prediction" %in% names(sub_copykat_res)) {
    pred <- as.data.frame(sub_copykat_res$prediction)
    pred$study_id <- current_id
    pred_list[[current_id]] <- pred
    qsave(pred, file.path(copykat_qs_dir, paste0("copykat_pred_", current_id, ".qs")))
    write.csv(pred, file.path(tables_dir, paste0("copykat_pred_", current_id, ".csv")), row.names = FALSE)
  } else {
    message("No valid CopyKAT prediction returned for ", current_id)
  }

  rm(sub_counts, sub_copykat_res)
  gc()
}

if (length(pred_list) == 0) {
  stop("No CopyKAT predictions were generated")
}

message("Step 4: merging CopyKAT predictions")
pred_merged <- bind_rows(pred_list)
qsave(pred_merged, file.path(copykat_qs_dir, "LUAD_copykat_merged_pred.qs"))
write.csv(pred_merged, file.path(tables_dir, "LUAD_copykat_merged_pred.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Add CopyKAT result to metadata and save malignant cells
# ------------------------------------------------------------------------------
message("Step 5: adding CopyKAT prediction to Seurat metadata")
pred_merged <- pred_merged %>%
  mutate(copykat.pred = as.character(copykat.pred))

copykat_status <- pred_merged$copykat.pred
names(copykat_status) <- pred_merged$cell.names
scRNA$copykat_pred <- copykat_status[colnames(scRNA)]
scRNA$copykat_pred[is.na(scRNA$copykat_pred)] <- "not_predicted"

copykat_summary <- scRNA@meta.data %>%
  dplyr::count(study_id, source_id, source_dataset, sample_group, copykat_pred, name = "cell_number")
write.csv(copykat_summary, file.path(tables_dir, "LUAD_copykat_summary_by_sample.csv"), row.names = FALSE)
qsave(scRNA, file.path(seurat_qs_dir, "LUAD_all_samples_RNA_harmony_copykat.qs"))

malig_barcodes <- pred_merged$cell.names[pred_merged$copykat.pred == "aneuploid"]
malig_barcodes <- intersect(malig_barcodes, colnames(scRNA))

malig_df <- scRNA@meta.data[malig_barcodes, , drop = FALSE] %>%
  tibble::rownames_to_column("Barcode")
write.table(
  malig_df["Barcode"],
  file = file.path(tables_dir, "LUAD_Malignant_cells.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.csv(malig_df, file.path(tables_dir, "LUAD_Malignant_cells_with_metadata.csv"), row.names = FALSE)

malig_counts <- malig_df %>%
  dplyr::count(study_id, source_id, source_dataset, sample_group, name = "number")

clinical <- fread(clinical_file, data.table = FALSE, na.strings = c("", "NA"))
clinical <- clinical %>%
  mutate(
    study_id = .data[["Patients ID in this study"]],
    source_id = .data[["Patients ID in GSE131907 / GSE123904"]]
  )

paint_malignant <- clinical %>%
  left_join(malig_counts, by = c("study_id", "source_id")) %>%
  mutate(number = if_else(is.na(number), 0L, as.integer(number))) %>%
  select(study_id, source_id, source_dataset, sample_group, number, everything())

write.table(
  paint_malignant,
  file = file.path(tables_dir, "LUAD_Paint_Malignant_cells.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.csv(paint_malignant, file.path(tables_dir, "LUAD_Paint_Malignant_cells.csv"), row.names = FALSE)

p_copykat <- DimPlot(scRNA, reduction = "umap", group.by = "copykat_pred", raster = TRUE)
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_copykat_pred.pdf"), p_copykat, width = 8, height = 6)
ggsave(file.path(plots_dir, "LUAD_all_samples_umap_by_copykat_pred.png"), p_copykat, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 5. Re-analyze malignant subset
# ------------------------------------------------------------------------------
message("Step 6: re-analyzing malignant subset")
if (length(malig_barcodes) < 50) {
  warning("Too few malignant cells for re-analysis: ", length(malig_barcodes))
} else {
  malignant <- subset(scRNA, cells = malig_barcodes)
  malignant <- JoinLayers(malignant)
  DefaultAssay(malignant) <- "RNA"

  malignant <- NormalizeData(malignant, normalization.method = "LogNormalize", scale.factor = 10000)
  malignant <- FindVariableFeatures(malignant, selection.method = "vst", nfeatures = 2000)
  malignant <- ScaleData(malignant, vars.to.regress = "percent.mt", verbose = FALSE)
  malignant <- RunPCA(malignant, npcs = 50, verbose = FALSE)
  malignant <- RunHarmony(
    malignant,
    group.by.vars = batch_var,
    assay.use = "RNA",
    max.iter.harmony = 20
  )

  malig_pcs <- 1:20
  malignant <- RunTSNE(malignant, reduction = "harmony", dims = malig_pcs)
  malignant <- RunUMAP(malignant, reduction = "harmony", dims = malig_pcs)
  malignant <- FindNeighbors(malignant, reduction = "harmony", dims = malig_pcs)
  malignant <- FindClusters(malignant, resolution = 1.0)

  qsave(malignant, file.path(seurat_qs_dir, "LUAD_Malignant_RNA_harmony.qs"))

  p_malig_cluster <- DimPlot(malignant, reduction = "umap", label = TRUE, raster = TRUE) + NoLegend()
  p_malig_sample <- DimPlot(malignant, reduction = "umap", group.by = "study_id", raster = TRUE)
  ggsave(file.path(plots_dir, "LUAD_malignant_umap_by_cluster.pdf"), p_malig_cluster, width = 8, height = 6)
  ggsave(file.path(plots_dir, "LUAD_malignant_umap_by_study_id.pdf"), p_malig_sample, width = 8, height = 6)
  ggsave(file.path(plots_dir, "LUAD_malignant_umap_by_cluster.png"), p_malig_cluster, width = 8, height = 6, dpi = 300)
  ggsave(file.path(plots_dir, "LUAD_malignant_umap_by_study_id.png"), p_malig_sample, width = 8, height = 6, dpi = 300)
}

message("Done. CopyKAT outputs are saved under: ", WORK_DIR)
