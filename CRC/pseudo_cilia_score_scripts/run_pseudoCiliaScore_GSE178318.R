suppressPackageStartupMessages({
  library(Seurat)
  library(qs)
})

source("/home/zerui/code/CRC/pseudo_cilia_score_scripts/pseudoCiliaScore_run_helpers.R")

DATASET <- "GSE178318"
WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE178318"
OUT_ROOT <- "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE178318"
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

setwd(WORK_DIR)
obj <- qread("Malignant.qs")
obj$type <- sub(".*_", "", colnames(obj))
obj <- subset(obj, subset = type %in% c("CRC", "LM"))
obj$group <- ifelse(obj$type == "LM", "Metastasis", "Primary")

cli <- read.csv("GSE178318_Cli.csv", header = TRUE, check.names = FALSE)
stage_df <- unique(data.frame(
  sample_id = cli[["Patient ID"]],
  Stage = cli[["Stage AJCC"]],
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
