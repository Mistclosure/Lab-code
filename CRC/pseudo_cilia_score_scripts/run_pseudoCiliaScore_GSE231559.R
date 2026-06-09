suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

source("/home/zerui/code/CRC/pseudo_cilia_score_scripts/pseudoCiliaScore_run_helpers.R")

DATASET <- "GSE231559"
WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE231559"
OUT_ROOT <- "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE231559"
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

setwd(WORK_DIR)
obj <- qread("Malignant_RNA_assay.qs")
cli <- read.csv("GSE231559_Cli.csv", header = TRUE, check.names = FALSE)

meta <- obj@meta.data
meta$cell_id_tmp <- rownames(meta)
meta_merged <- merge(meta, cli, by = "Patient", all.x = TRUE, sort = FALSE)
rownames(meta_merged) <- meta_merged$cell_id_tmp
meta_merged$cell_id_tmp <- NULL
obj@meta.data <- meta_merged[colnames(obj), , drop = FALSE]
obj$group <- as.character(obj@meta.data[["metastasis"]])
obj$Stage <- as.character(obj@meta.data[["Stage"]])

stage_df <- unique(data.frame(
  sample_id = obj@meta.data[["orig.ident"]],
  Stage = obj@meta.data[["Stage"]],
  stringsAsFactors = FALSE
))

run_pseudo_cilia_for_gene_sets(
  obj = obj,
  dataset = DATASET,
  out_root = OUT_ROOT,
  sample_col = "orig.ident",
  group_col = "group",
  stage_df = stage_df,
  stage_colname = "Stage",
  cells_per_pseudo = 80,
  min_cells_per_pseudo = 30
)

message("Finished ", DATASET, " pseudoCiliaScore runs.")
