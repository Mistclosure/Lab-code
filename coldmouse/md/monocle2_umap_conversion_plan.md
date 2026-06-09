# Pseudotime Trajectory Analysis_standard.r 转 monocle2 + UMAP 展示方案

目标脚本：

`/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_standard.r`

目标：

将当前 monocle3 版本改为 monocle2 版本。保留原始 PBMC monocyte 与 Aorta macrophage 的筛选、合并、分组过滤和 metadata 标注逻辑。拟时序计算使用 monocle2，最终图使用 Seurat UMAP 坐标展示 monocle2 的 `Pseudotime`、`State`、`Source` 和 `Group`，避免输出 monocle2 默认的 component 图。

重要说明：

monocle2 标准轨迹推断依赖 `reduceDimension(method = "DDRTree")` 和 `orderCells()`。monocle2 不能像 monocle3 一样直接用 UMAP 学习 principal graph。因此不要简单把 `reduceDimension(method = "DDRTree")` 改成 UMAP。正确改法是：

1. monocle2 内部仍用 DDRTree 完成拟时序排序。
2. 不使用 `plot_cell_trajectory()` 作为最终主图。
3. 将 monocle2 计算得到的 `Pseudotime` 和 `State` 合并到 Seurat metadata。
4. 使用 Seurat 已有 UMAP 坐标，通过 `ggplot2` 绘制 UMAP 上的 monocle2 结果。

## 1. 包依赖修改

删除：

```r
library(monocle3)
library(harmony)
```

改为：

```r
library(monocle)
library(Matrix)
library(Biobase)
```

保留：

```r
library(qs)
library(Seurat)
library(ggplot2)
library(dplyr)
```

如果运行时 monocle2 与新版 `dplyr` 或 `igraph` 有兼容问题，参考：

`/home/zerui/code/coldmouse_introns/scripts/07_monocle2_trajectory_introns.R`

只复制其中 monocle2 compatibility patch，避免重构整个脚本。

## 2. 保留原数据读取和细胞筛选逻辑

以下逻辑保持不变：

```r
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

pbmc <- qread('pbmc_corrected.qs')
aorta <- qread('aorta_corrected.qs')

pbmc_sub <- subset(pbmc, new_clusters %in% c('3', '12'))
pbmc_sub$cell_type <- ifelse(
  pbmc_sub$new_clusters == '12',
  "Monocyte",
  "Dendritic cells"
)
pbmc_sub <- subset(pbmc_sub, new_clusters == '12')

aorta_target <- subset(aorta, seurat_clusters %in% c('1','6','10','14'))
aorta_target$cell_type <- aorta_target$SingleR.labels

merged_obj <- merge(
  pbmc_sub,
  y = aorta_target,
  add.cell.ids = c("PBMC", "Aorta")
)

merged_obj <- subset(merged_obj, Group == 'TN_30C', invert = TRUE)
merged_obj$celllabels <- sub(".*:", "", merged_obj$celltype)
merged_obj$Source <- merged_obj$celllabels
merged_obj <- JoinLayers(merged_obj)

print(table(merged_obj$Source, merged_obj$cell_type))
```

## 3. 确保 Seurat UMAP 可用

在构建 monocle2 对象前检查 `merged_obj` 是否已有 UMAP：

```r
if (!"umap" %in% names(merged_obj@reductions)) {
  message("No UMAP reduction found in merged_obj; running Seurat UMAP.")
  DefaultAssay(merged_obj) <- "RNA"
  merged_obj <- NormalizeData(merged_obj)
  merged_obj <- FindVariableFeatures(merged_obj, nfeatures = 2000)
  merged_obj <- ScaleData(merged_obj)
  merged_obj <- RunPCA(merged_obj, npcs = 30)
  merged_obj <- RunUMAP(merged_obj, dims = 1:30)
}

umap_coords <- Embeddings(merged_obj, "umap")
```

说明：

如果原始 `pbmc_corrected.qs` 和 `aorta_corrected.qs` 合并后 UMAP 不可用，就现场用 Seurat 重新计算 UMAP。这个 UMAP 只用于最终展示，不用于 monocle2 的拟时序排序。

## 4. 构建 monocle2 CellDataSet

替换原 monocle3 的：

```r
new_cell_data_set()
estimate_size_factors()
preprocess_cds()
align_cds()
reduce_dimension()
cluster_cells()
learn_graph()
order_cells()
```

改为 monocle2：

```r
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
```

## 5. 设置 ordering genes

优先使用 Seurat 的 variable features。若数量不足，则回退到 monocle2 检测到的 expressed genes：

```r
expressed_genes <- rownames(subset(fData(cds), num_cells_expressed >= 10))

ordering_genes <- VariableFeatures(merged_obj)
ordering_genes <- ordering_genes[ordering_genes %in% expressed_genes]

if (length(ordering_genes) < 100) {
  warning("VariableFeatures < 100; fallback to expressed genes.")
  ordering_genes <- expressed_genes
}

cds <- setOrderingFilter(cds, ordering_genes)
```

## 6. monocle2 拟时序排序

注意：这一步仍需 DDRTree，因为这是 monocle2 的核心排序机制。最终展示会改用 UMAP，不输出 DDRTree component 主图。

```r
cds <- reduceDimension(
  cds,
  max_components = 2,
  method = "DDRTree"
)

cds <- orderCells(cds)
```

## 7. 使用 Ly6c2 高表达细胞确定 root_state

替换 monocle3 的 `root_pr_nodes` 逻辑，改成 monocle2 的 `root_state`：

```r
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
cds <- orderCells(cds, root_state = root_state)
```

## 8. 将 monocle2 结果合并回 Seurat metadata

```r
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
```

## 9. 使用 UMAP 绘制 monocle2 结果

不要再用 `plot_cell_trajectory()` 作为主输出图。改用 Seurat UMAP 坐标：

```r
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
```

## 10. 保存输出

输出文件名必须使用 `monocle2_umap`，避免和 monocle3 或 DDRTree component 图混淆：

```r
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
```

## 给 Claude Code 的执行提示词

```text
请修改以下脚本，不要运行完整分析：

/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_standard.r

目标：
把当前 monocle3 workflow 改成 monocle2 workflow，并且最终结果图使用 Seurat UMAP 坐标展示 monocle2 的 Pseudotime 和 State，不要输出 monocle2 默认 component 图作为主结果。

具体要求：
1. 保留原始数据路径、pbmc_corrected.qs、aorta_corrected.qs 读取逻辑。
2. 保留 PBMC new_clusters == 12 的 Monocyte 筛选。
3. 保留 Aorta seurat_clusters in c('1','6','10','14') 的 Macrophage 相关筛选。
4. 保留 merge、Group != TN_30C、celllabels、Source、cell_type 的原逻辑。
5. 删除 monocle3 和 harmony 依赖，改用 monocle、Matrix、Biobase。
6. 用 LayerData(merged_obj, assay = "RNA", layer = "counts") 获取 counts。
7. 用 AnnotatedDataFrame 构建 monocle2 的 phenoData 和 featureData。
8. 用 newCellDataSet()、estimateSizeFactors()、estimateDispersions()、detectGenes() 构建 monocle2 对象。
9. ordering genes 优先使用 VariableFeatures(merged_obj)，如果少于 100 个则回退到 num_cells_expressed >= 10 的 genes。
10. monocle2 拟时序排序内部仍使用 reduceDimension(method = "DDRTree") 和 orderCells()，因为 monocle2 不支持直接用 UMAP 学习轨迹。不要把 method 改成 UMAP。
11. 用 Ly6c2 高表达细胞所在最多的 State 作为 root_state，并再次 orderCells(cds, root_state = root_state)。
12. 不使用 plot_cell_trajectory() 作为主图输出。
13. 将 pData(cds)$Pseudotime 和 pData(cds)$State 合并回 merged_obj metadata。
14. 如果 merged_obj 没有 umap reduction，则用 Seurat 重新 NormalizeData、FindVariableFeatures、ScaleData、RunPCA、RunUMAP 生成 UMAP。
15. 使用 Embeddings(merged_obj, "umap") 和 ggplot2 绘制：
    - monocle2_umap_source
    - monocle2_umap_pseudotime
    - monocle2_umap_state
    - monocle2_umap_split
16. 保存 csv、qs、rds、png、pdf。输出文件名必须包含 monocle2_umap，避免与 monocle3 或 DDRTree component 图混淆。
17. 修改完成后只做 R 语法检查，不运行完整分析。
18. 如果 monocle2 和新版 dplyr/igraph 出现兼容问题，参考 /home/zerui/code/coldmouse_introns/scripts/07_monocle2_trajectory_introns.R 中的 compatibility patch，最小范围复制。
```
