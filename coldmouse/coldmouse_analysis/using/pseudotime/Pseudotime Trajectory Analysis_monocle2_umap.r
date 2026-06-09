# ==============================================================================
# Pseudotime Trajectory Analysis_monocle2_umap.r
# Monocle 2 workflow with Seurat UMAP visualization
# PBMC Monocyte + Aorta Macrophage combined pseudotime analysis
# ==============================================================================
# Monocle2 内部仍用 DDRTree 完成拟时序排序。
# 最终图使用 Seurat UMAP 坐标展示 monocle2 的 Pseudotime、State、Source 和 Group。
# 不输出 monocle2 默认的 component 图。
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

library(qs)
library(Seurat)
library(harmony)
library(ggplot2)
library(dplyr)
library(Matrix)
library(Biobase)

# ------------------------------------------------------------------------------
# monocle2 compatibility patches for current dplyr / igraph
# ------------------------------------------------------------------------------
# Patch dplyr namespace: monocle calls dplyr::group_by_() etc. directly
for (fn_name in c("group_by_", "summarise_", "mutate_", "filter_", "arrange_",
                  "select_", "slice_", "distinct_")) {
  base_fn <- gsub("_$", "", fn_name)
  fn_def <- eval(parse(text = paste0("function(.data, ...) { dplyr::", base_fn, "(.data, ...) }")))
  assignInNamespace(fn_name, fn_def, ns = "dplyr")
  if (!exists(fn_name, envir = globalenv())) {
    assign(fn_name, fn_def, envir = globalenv())
  }
}

library(monocle)

if (requireNamespace("igraph", quietly = TRUE)) {
  orig_dfs <- igraph::dfs
  patched_dfs <- function(graph, root, neimode = NULL, mode = NULL, ...) {
    if (!is.null(neimode) && is.null(mode)) mode <- neimode
    orig_dfs(graph, root, mode = mode, ...)
  }
  environment(patched_dfs) <- environment(orig_dfs)
  assignInNamespace("dfs", patched_dfs, ns = "igraph")
}

patch_monocle_nei <- function() {
  ns <- getNamespace("monocle")
  for (fn_name in ls(ns)) {
    fn <- get(fn_name, envir = ns)
    if (!is.function(fn)) next
    body_txt <- paste(deparse(body(fn)), collapse = "\n")
    if (!grepl("nei(", body_txt, fixed = TRUE)) next
    patched_txt <- gsub("(?<![.])nei\\(", ".nei(", body_txt, perl = TRUE)
    patched_body <- tryCatch(parse(text = patched_txt)[[1]], error = function(e) NULL)
    if (is.null(patched_body)) next
    body(fn) <- patched_body
    assignInNamespace(fn_name, fn, ns = "monocle")
  }
}
patch_monocle_nei()

# ==============================================================================
# 1. 读取数据
# ==============================================================================
print("Loading data...")
pbmc <- qread('pbmc_recorrected.qs')
aorta <- qread('aorta_corrected.qs')

# ==============================================================================
# 2. 提取与标注
# ==============================================================================
print("Subsetting PBMC Monocytes by cell_annotation...")
if (!"cell_annotation" %in% colnames(pbmc@meta.data)) {
  stop("PBMC object does not contain metadata column: cell_annotation")
}
pbmc_sub <- subset(pbmc, cell_annotation == "Monocytes")
pbmc_sub$cell_type <- "Monocytes"
pbmc_sub$Source <- "PBMC_Monocytes"
pbmc_sub$batch_source <- "PBMC"

print("Subsetting Aorta Macrophage/Monocyte cells by cell_annotation...")
if (!"cell_annotation" %in% colnames(aorta@meta.data)) {
  stop("Aorta object does not contain metadata column: cell_annotation")
}
aorta_target <- subset(aorta, grepl("Macrophage|Monocyte", cell_annotation, ignore.case = TRUE))
aorta_target$cell_type <- aorta_target$cell_annotation
aorta_target$Source <- "Aorta_Macrophage_Monocyte"
aorta_target$batch_source <- "Aorta"

# ==============================================================================
# 3. 合并与预处理
# ==============================================================================
print("Merging objects...")
merged_obj <- merge(pbmc_sub, y = aorta_target,
                    add.cell.ids = c("PBMC", "Aorta"))
merged_obj <- subset(merged_obj, Group == 'TN_30C', invert = TRUE)
merged_obj$celllabels <- merged_obj$cell_annotation
merged_obj$Source <- ifelse(grepl("^PBMC_", colnames(merged_obj)),
                            "PBMC_Monocytes", "Aorta_Macrophage_Monocyte")
merged_obj$batch_source <- ifelse(grepl("^PBMC_", colnames(merged_obj)),
                                  "PBMC", "Aorta")
merged_obj$cell_type <- ifelse(grepl("^PBMC_", colnames(merged_obj)),
                               "Monocytes", merged_obj$cell_annotation)
merged_obj <- JoinLayers(merged_obj)

print(table(merged_obj$Source, merged_obj$cell_type))
print(table(merged_obj$batch_source, merged_obj$Group))

# ==============================================================================
# 4. Harmony 批次校正 + Seurat UMAP
# ==============================================================================
print("Running Harmony integration to correct batch effects (PBMC vs Aorta)...")
DefaultAssay(merged_obj) <- "RNA"
merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj, npcs = 30)

merged_obj <- RunHarmony(merged_obj, group.by.vars = "batch_source", reduction.use = "pca",
                         dims.use = 1:30)

merged_obj <- RunUMAP(merged_obj, reduction = "harmony", dims = 1:30)

# ==============================================================================
# 5. 构建 monocle2 CellDataSet
# ==============================================================================
print("Converting to Monocle 2 CellDataSet...")
expression_matrix <- LayerData(merged_obj, assay = "RNA", layer = "counts")

cell_metadata <- merged_obj@meta.data
cell_metadata <- cell_metadata[colnames(expression_matrix), , drop = FALSE]

gene_metadata <- data.frame(
  gene_short_name = rownames(expression_matrix),
  row.names = rownames(expression_matrix)
)

pd <- new("AnnotatedDataFrame", data = cell_metadata)
fd <- new("AnnotatedDataFrame", data = gene_metadata)

cds <- newCellDataSet(
  as(expression_matrix, "sparseMatrix"),
  phenoData = pd,
  featureData = fd,
  expressionFamily = negbinomial.size()
)

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)

# ==============================================================================
# 6. 设置 ordering genes
# ==============================================================================
print("Setting ordering genes...")
expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

ordering_genes <- VariableFeatures(merged_obj)
ordering_genes <- ordering_genes[ordering_genes %in% expressed_genes]

if (length(ordering_genes) < 100) {
  warning("VariableFeatures < 100; fallback to expressed genes.")
  ordering_genes <- expressed_genes
}

cds <- setOrderingFilter(cds, ordering_genes)
print(paste0("Number of ordering genes: ", length(ordering_genes)))

# ==============================================================================
# 7. monocle2 拟时序排序 (DDRTree)
# ==============================================================================
print("Running reduceDimension with DDRTree...")
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree"
)

print("Ordering cells...")
cds <- orderCells(cds)

# ==============================================================================
# 8. 使用 Ly6c2 高表达细胞确定 root_state
# ==============================================================================
print("Determining root state by Ly6c2 expression...")
get_root_state_by_marker <- function(cds, marker_gene = "Ly6c2", q = 0.95) {
  if (!marker_gene %in% rownames(cds)) {
    warning(marker_gene, " not found; using largest State as root.")
    return(as.numeric(names(sort(table(pData(cds)$State), decreasing = TRUE))[1]))
  }

  gene_expr <- exprs(cds)[marker_gene, ]
  cutoff <- quantile(gene_expr, q, na.rm = TRUE)
  early_cells <- names(gene_expr[gene_expr >= cutoff])

  if (length(early_cells) == 0) {
    warning("No high-expression cells found for ", marker_gene, "; using largest State as root.")
    return(as.numeric(names(sort(table(pData(cds)$State), decreasing = TRUE))[1]))
  }

  root_state <- names(sort(table(pData(cds)[early_cells, "State"]), decreasing = TRUE))[1]
  as.numeric(root_state)
}

root_state <- get_root_state_by_marker(cds, marker_gene = "Ly6c2")
print(paste0("Root state: ", root_state))
cds <- orderCells(cds, root_state = root_state)

# ==============================================================================
# 9. 将 monocle2 结果合并回 Seurat metadata
# ==============================================================================
print("Merging monocle2 results back to Seurat object...")
pseudotime_meta <- data.frame(
  Cell = rownames(pData(cds)),
  Pseudotime = pData(cds)$Pseudotime,
  State = pData(cds)$State,
  Source = pData(cds)$Source,
  Group = pData(cds)$Group,
  cell_type = pData(cds)$cell_type,
  stringsAsFactors = FALSE
)

merged_obj$monocle2_pseudotime <- NA_real_
merged_obj$monocle2_state <- NA

merged_obj$monocle2_pseudotime[pseudotime_meta$Cell] <- pseudotime_meta$Pseudotime
merged_obj$monocle2_state[pseudotime_meta$Cell] <- pseudotime_meta$State

# ==============================================================================
# 10. 使用 UMAP 绘制 monocle2 结果
# ==============================================================================
print("Generating UMAP plots...")
umap_coords <- Embeddings(merged_obj, "umap")

umap_df <- data.frame(
  Cell = rownames(umap_coords),
  UMAP_1 = umap_coords[, 1],
  UMAP_2 = umap_coords[, 2],
  Source = merged_obj$Source,
  Group = merged_obj$Group,
  cell_type = merged_obj$cell_type,
  monocle2_pseudotime = merged_obj$monocle2_pseudotime,
  monocle2_state = merged_obj$monocle2_state,
  stringsAsFactors = FALSE
)

plot_umap_continuous <- function(df, color_by, title) {
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]])) +
    geom_point(size = 0.35, alpha = 0.85) +
    scale_color_gradientn(
      colors = c("#2c7bb6", "#ffffbf", "#d7191c"),
      na.value = "grey80"
    ) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2", color = color_by) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_umap_discrete <- function(df, color_by, title) {
  df[[color_by]] <- factor(df[[color_by]])
  ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]])) +
    geom_point(size = 0.35, alpha = 0.85) +
    labs(title = title, x = "UMAP_1", y = "UMAP_2", color = color_by) +
    theme_classic(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5))
}

plot_source_umap <- plot_umap_discrete(
  umap_df,
  "Source",
  "Monocle2 Result on Seurat UMAP: Source"
)

plot_pseudotime_umap <- plot_umap_continuous(
  umap_df,
  "monocle2_pseudotime",
  "Monocle2 Pseudotime on Seurat UMAP"
)

plot_state_umap <- plot_umap_discrete(
  umap_df,
  "monocle2_state",
  "Monocle2 State on Seurat UMAP"
)

plot_pseudotime_split_umap <- plot_pseudotime_umap +
  facet_wrap(~Group) +
  ggtitle("Monocle2 Pseudotime on Seurat UMAP: Cold vs RT")

# ==============================================================================
# 11. 保存输出
# ==============================================================================
print("Saving outputs...")
if (!dir.exists("files")) dir.create("files")
if (!dir.exists("pictures")) dir.create("pictures")

qsave(cds, "files/script_06_monocle2_standard.qs")
saveRDS(cds, "files/script_06_monocle2_standard.rds")
qsave(merged_obj, "files/script_06_monocle2_umap_seurat_with_pseudotime.qs")

write.csv(
  pseudotime_meta,
  "files/script_06_monocle2_pseudotime_metadata.csv",
  row.names = FALSE
)

write.csv(
  umap_df,
  "files/script_06_monocle2_umap_metadata.csv",
  row.names = FALSE
)

ggsave("pictures/monocle2_umap_source.png", plot_source_umap, width = 8, height = 5, dpi = 300)
ggsave("pictures/monocle2_umap_pseudotime.png", plot_pseudotime_umap, width = 6, height = 5, dpi = 300)
ggsave("pictures/monocle2_umap_state.png", plot_state_umap, width = 6, height = 5, dpi = 300)
ggsave("pictures/monocle2_umap_split.png", plot_pseudotime_split_umap, width = 10, height = 5, dpi = 300)

ggsave("pictures/monocle2_umap_source.pdf", plot_source_umap, width = 8, height = 5)
ggsave("pictures/monocle2_umap_pseudotime.pdf", plot_pseudotime_umap, width = 6, height = 5)
ggsave("pictures/monocle2_umap_state.pdf", plot_state_umap, width = 6, height = 5)
ggsave("pictures/monocle2_umap_split.pdf", plot_pseudotime_split_umap, width = 10, height = 5)

print("Monocle2 UMAP analysis completed successfully!")
