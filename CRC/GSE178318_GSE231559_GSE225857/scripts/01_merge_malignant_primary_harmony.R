suppressPackageStartupMessages({
  library(data.table)
  library(Seurat)
  library(harmony)
  library(qs)
  library(dplyr)
  library(ggplot2)
})

set.seed(1234)

base_dir <- "/mnt/disk1/qiuzerui/downloads/CRC"
project_name <- "GSE178318_GSE231559_GSE225857"
out_dir <- file.path(base_dir, project_name)
files_dir <- file.path(out_dir, "files")
plots_dir <- file.path(out_dir, "plots")
qs_dir <- file.path(out_dir, "qs")
logs_dir <- file.path(out_dir, "logs")

dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(qs_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(logs_dir, "01_merge_malignant_primary_harmony.log")
sink(log_file, split = TRUE)
on.exit(sink(), add = TRUE)

message("Output directory: ", out_dir)

read_malignant <- function(path, dataset) {
  if (!file.exists(path)) {
    stop("Input qs not found: ", path)
  }
  obj <- qread(path)
  DefaultAssay(obj) <- "RNA"
  if ("SCT" %in% Assays(obj)) {
    obj[["SCT"]] <- NULL
  }
  obj <- JoinLayers(obj)
  obj$source_dataset <- dataset
  obj
}

add_selection_metadata <- function(obj, clinical_info, sample_col = "orig.ident") {
  md <- obj@meta.data
  key <- md[[sample_col]]
  rownames(clinical_info) <- clinical_info$source_sample_raw

  obj$source_sample_raw <- key
  obj$patient_id_raw <- clinical_info[key, "patient_id_raw"]
  obj$source_sample <- clinical_info[key, "source_sample"]
  obj$patient_id <- clinical_info[key, "patient_id"]
  obj$metastasis_status <- clinical_info[key, "metastasis_status"]
  obj$primary_label <- clinical_info[key, "primary_label"]
  obj$lesion_site <- clinical_info[key, "lesion_site"]
  obj$clinical_stage <- clinical_info[key, "clinical_stage"]
  obj$selection_rule <- clinical_info[key, "selection_rule"]
  obj
}

rebuild_from_rna_counts <- function(obj, project) {
  DefaultAssay(obj) <- "RNA"
  counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
  CreateSeuratObject(
    counts = counts,
    meta.data = obj@meta.data,
    project = project
  )
}

message("Step 1: reading malignant Seurat objects")
gse178318 <- read_malignant(
  file.path(base_dir, "GSE178318", "Malignant.qs"),
  "GSE178318"
)
gse231559 <- read_malignant(
  file.path(base_dir, "GSE231559", "Malignant_GSE231559.qs"),
  "GSE231559"
)
gse225857 <- read_malignant(
  file.path(base_dir, "GSE225857", "Malignant.qs"),
  "GSE225857"
)

message("Input cells:")
print(table(gse178318$orig.ident, useNA = "ifany"))
print(table(gse231559$orig.ident, useNA = "ifany"))
print(table(gse225857$orig.ident, useNA = "ifany"))

message("Step 2: reading clinical tables and selecting primary malignant cells")
cli178318 <- fread(file.path(base_dir, "GSE178318", "GSE178318_Cli.csv"), data.table = FALSE)
cli231559 <- fread(file.path(base_dir, "GSE231559", "GSE231559_Cli.csv"), data.table = FALSE)
cli225857 <- fread(file.path(base_dir, "GSE225857", "GSE225857_Cli.csv"), data.table = FALSE)

clin178318 <- cli178318 %>%
  transmute(
    source_dataset = "GSE178318",
    source_sample_raw = .data[["Patient ID"]],
    source_sample = paste0("GSE178318_", .data[["Patient ID"]]),
    patient_id_raw = .data[["Patient ID"]],
    patient_id = paste0("GSE178318_", .data[["Patient ID"]]),
    metastasis_status = "metastatic_patient",
    primary_label = "metastatic_patient_primary",
    lesion_site = "primary",
    clinical_stage = .data[["Stage AJCC"]],
    selection_rule = "All malignant cells in GSE178318 are CRC primary tumors from patients with liver metastasis."
  )

primary231559 <- cli231559 %>%
  filter(.data[["metastasis"]] == "Primary")
clin231559 <- primary231559 %>%
  transmute(
    source_dataset = "GSE231559",
    source_sample_raw = paste0(.data[["Patient"]], "T"),
    source_sample = paste0("GSE231559_", .data[["Patient"]], "T"),
    patient_id_raw = .data[["Patient"]],
    patient_id = paste0("GSE231559_", .data[["Patient"]]),
    metastasis_status = "non_metastatic_patient",
    primary_label = "non_metastatic_patient_primary",
    lesion_site = "primary",
    clinical_stage = .data[["Stage"]],
    selection_rule = "Primary samples in GSE231559 clinical table (metastasis == Primary)."
  )

primary225857 <- cli225857 %>%
  filter(.data[["Colorectal cancer (CC)"]] == "Yes")
clin225857 <- primary225857 %>%
  transmute(
    source_dataset = "GSE225857",
    source_sample_raw = .data[["Patient ID"]],
    source_sample = paste0("GSE225857_", .data[["Patient ID"]]),
    patient_id_raw = .data[["Patient ID"]],
    patient_id = paste0("GSE225857_", .data[["Patient ID"]]),
    metastasis_status = "non_metastatic_patient",
    primary_label = "non_metastatic_patient_primary",
    lesion_site = "primary",
    clinical_stage = NA_character_,
    selection_rule = "Primary colorectal cancer samples marked CC=Yes; LM samples are excluded."
  )

clinical_selected <- bind_rows(clin178318, clin231559, clin225857)
write.csv(
  clinical_selected %>% select(patient_id, patient_id_raw, metastasis_status, primary_label),
  file.path(files_dir, "clinical_info_minimal.csv"),
  row.names = FALSE
)
write.csv(
  clinical_selected,
  file.path(files_dir, "clinical_info_selected_samples.csv"),
  row.names = FALSE
)

selected178318 <- unique(clin178318$source_sample_raw)
selected231559 <- unique(clin231559$source_sample_raw)
selected225857 <- unique(clin225857$source_sample_raw)

gse178318 <- subset(gse178318, subset = orig.ident %in% selected178318)
gse231559 <- subset(gse231559, subset = orig.ident %in% selected231559)
gse225857 <- subset(gse225857, subset = orig.ident %in% selected225857)

if (ncol(gse178318) == 0 || ncol(gse231559) == 0 || ncol(gse225857) == 0) {
  stop(
    "At least one selected object has zero cells: ",
    "GSE178318=", ncol(gse178318), ", ",
    "GSE231559=", ncol(gse231559), ", ",
    "GSE225857=", ncol(gse225857)
  )
}

gse178318 <- add_selection_metadata(gse178318, clin178318)
gse231559 <- add_selection_metadata(gse231559, clin231559)
gse225857 <- add_selection_metadata(gse225857, clin225857)

gse178318 <- rebuild_from_rna_counts(gse178318, "GSE178318")
gse231559 <- rebuild_from_rna_counts(gse231559, "GSE231559")
gse225857 <- rebuild_from_rna_counts(gse225857, "GSE225857")

gse178318 <- RenameCells(gse178318, add.cell.id = "GSE178318")
gse231559 <- RenameCells(gse231559, add.cell.id = "GSE231559")
gse225857 <- RenameCells(gse225857, add.cell.id = "GSE225857")

message("Selected cells:")
selected_summary <- bind_rows(
  gse178318@meta.data %>% count(source_dataset, source_sample, source_sample_raw, patient_id, patient_id_raw, metastasis_status, primary_label, name = "cell_number"),
  gse231559@meta.data %>% count(source_dataset, source_sample, source_sample_raw, patient_id, patient_id_raw, metastasis_status, primary_label, name = "cell_number"),
  gse225857@meta.data %>% count(source_dataset, source_sample, source_sample_raw, patient_id, patient_id_raw, metastasis_status, primary_label, name = "cell_number")
)
print(selected_summary)
write.csv(selected_summary, file.path(files_dir, "selected_cell_counts_by_sample.csv"), row.names = FALSE)

message("Step 3: merging selected malignant primary cells")
combined <- merge(
  gse178318,
  y = list(gse231559, gse225857),
  project = project_name,
  merge.data = FALSE
)
combined <- JoinLayers(combined)
DefaultAssay(combined) <- "RNA"

mt_pattern <- if (any(grepl("^MT-", rownames(combined)))) "^MT-" else "^mt-"
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = mt_pattern)

message("Merged object: ", nrow(combined), " genes x ", ncol(combined), " cells")
print(table(combined$source_dataset, combined$primary_label, useNA = "ifany"))

message("Step 4: RNA workflow, Harmony integration, UMAP, and clustering")
combined <- NormalizeData(combined, normalization.method = "LogNormalize", scale.factor = 10000)
combined <- FindVariableFeatures(combined, selection.method = "vst", nfeatures = 2000)
combined <- ScaleData(combined, vars.to.regress = "percent.mt", verbose = FALSE)
combined <- RunPCA(combined, npcs = 50, verbose = FALSE)
combined <- RunHarmony(
  combined,
  group.by.vars = "source_dataset",
  assay.use = "RNA",
  max.iter.harmony = 20
)

pcs_use <- 1:20
combined <- RunTSNE(combined, reduction = "harmony", dims = pcs_use)
combined <- RunUMAP(combined, reduction = "harmony", dims = pcs_use)
combined <- FindNeighbors(combined, reduction = "harmony", dims = pcs_use)
combined <- FindClusters(combined, resolution = 1.0)

message("Step 5: saving qs, metadata, and plots")
qsave(combined, file.path(qs_dir, "GSE178318_GSE231559_GSE225857_malignant_primary_RNA_harmony.qs"))
write.csv(
  combined@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    select(cell_id, patient_id, patient_id_raw, metastasis_status, primary_label, source_dataset, source_sample, source_sample_raw, lesion_site, clinical_stage, seurat_clusters),
  file.path(files_dir, "merged_malignant_primary_metadata.csv"),
  row.names = FALSE
)

p_cluster <- DimPlot(combined, reduction = "umap", label = TRUE, raster = TRUE) + NoLegend()
p_dataset <- DimPlot(combined, reduction = "umap", group.by = "source_dataset", raster = TRUE)
p_sample <- DimPlot(combined, reduction = "umap", group.by = "source_sample", raster = TRUE)
p_label <- DimPlot(combined, reduction = "umap", group.by = "primary_label", raster = TRUE)
p_metastasis <- DimPlot(combined, reduction = "umap", group.by = "metastasis_status", raster = TRUE)

ggsave(file.path(plots_dir, "malignant_primary_umap_by_cluster.pdf"), p_cluster, width = 8, height = 6)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_dataset.pdf"), p_dataset, width = 8, height = 6)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_sample.pdf"), p_sample, width = 9, height = 6)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_primary_label.pdf"), p_label, width = 8, height = 6)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_metastasis_status.pdf"), p_metastasis, width = 8, height = 6)

ggsave(file.path(plots_dir, "malignant_primary_umap_by_cluster.png"), p_cluster, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_dataset.png"), p_dataset, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_sample.png"), p_sample, width = 9, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_primary_label.png"), p_label, width = 8, height = 6, dpi = 300)
ggsave(file.path(plots_dir, "malignant_primary_umap_by_metastasis_status.png"), p_metastasis, width = 8, height = 6, dpi = 300)

cluster_summary <- combined@meta.data %>%
  count(source_dataset, source_sample, source_sample_raw, patient_id, patient_id_raw, metastasis_status, primary_label, seurat_clusters, name = "cell_number")
write.csv(cluster_summary, file.path(files_dir, "cluster_cell_counts_by_sample.csv"), row.names = FALSE)

message("Done.")
