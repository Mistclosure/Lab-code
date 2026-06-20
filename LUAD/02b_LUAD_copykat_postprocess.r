# ==============================================================================
# LUAD CopyKAT post-processing
# Continue from existing CopyKAT predictions without re-running CopyKAT.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(harmony)
  library(qs)
  library(dplyr)
  library(ggplot2)
  library(tibble)
})

set.seed(1234)

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/LUAD"
setwd(WORK_DIR)

input_qs <- file.path(WORK_DIR, "qs", "Seurat", "LUAD_all_samples_RNA_harmony_before_copykat.qs")
pred_qs <- file.path(WORK_DIR, "qs", "copykat", "LUAD_copykat_merged_pred.qs")
clinical_file <- file.path(WORK_DIR, "LUAD_clinical_merged_from_image.csv")

tables_dir <- file.path(WORK_DIR, "files", "tables")
plots_dir <- file.path(WORK_DIR, "plots", "copykat")
seurat_qs_dir <- file.path(WORK_DIR, "qs", "Seurat")
logs_dir <- file.path(WORK_DIR, "logs")

dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(seurat_qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "02b_LUAD_copykat_postprocess.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("Step 1: reading integrated Seurat object and CopyKAT predictions")
if (!file.exists(input_qs)) {
  stop("Input Seurat qs not found: ", input_qs)
}
if (!file.exists(pred_qs)) {
  stop("Merged CopyKAT prediction qs not found: ", pred_qs)
}
if (!file.exists(clinical_file)) {
  stop("Clinical file not found: ", clinical_file)
}

scRNA <- qread(input_qs)
pred_merged <- qread(pred_qs)

required_pred_cols <- c("cell.names", "copykat.pred", "study_id")
missing_pred_cols <- setdiff(required_pred_cols, colnames(pred_merged))
if (length(missing_pred_cols) > 0) {
  stop("Missing CopyKAT prediction columns: ", paste(missing_pred_cols, collapse = ", "))
}

message("Seurat object: ", nrow(scRNA), " genes x ", ncol(scRNA), " cells")
message("CopyKAT predicted cells: ", nrow(pred_merged))
print(table(pred_merged$copykat.pred, useNA = "ifany"))

message("Step 2: adding CopyKAT prediction to Seurat metadata")
pred_merged <- pred_merged %>%
  mutate(copykat.pred = as.character(copykat.pred))

copykat_status <- pred_merged$copykat.pred
names(copykat_status) <- pred_merged$cell.names
scRNA$copykat_pred <- unname(copykat_status[colnames(scRNA)])
scRNA$copykat_pred[is.na(scRNA$copykat_pred)] <- "not_predicted"

copykat_summary <- scRNA@meta.data %>%
  dplyr::count(study_id, source_id, source_dataset, sample_group, copykat_pred, name = "cell_number")
write.csv(copykat_summary, file.path(tables_dir, "LUAD_copykat_summary_by_sample.csv"), row.names = FALSE)
qsave(scRNA, file.path(seurat_qs_dir, "LUAD_all_samples_RNA_harmony_copykat.qs"))

message("Step 3: exporting malignant-cell metadata and clinical count table")
malig_barcodes <- pred_merged$cell.names[pred_merged$copykat.pred == "aneuploid"]
malig_barcodes <- intersect(malig_barcodes, colnames(scRNA))
message("Malignant cells retained in Seurat object: ", length(malig_barcodes))

malig_df <- scRNA@meta.data[malig_barcodes, , drop = FALSE] %>%
  tibble::rownames_to_column("cell_barcode")
write.table(
  data.frame(Barcode = malig_df$cell_barcode),
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

message("Step 4: re-analyzing malignant subset")
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
    group.by.vars = "study_id",
    assay.use = "RNA",
    max_iter = 20
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

message("Done. CopyKAT post-processing outputs are saved under: ", WORK_DIR)
