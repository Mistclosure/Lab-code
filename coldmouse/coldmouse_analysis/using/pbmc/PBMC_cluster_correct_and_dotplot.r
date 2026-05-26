# 此脚本用于纠正分群后的出图，请在02之后使用
library(Seurat)
library(qs)
library(patchwork)
library(ggplot2)
library(dplyr)

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

sc_by_tissue = qread('sc_by_tissue_louvain.qs')
pbmc = sc_by_tissue[['PBMC']]

# 1. 剔除 cluster 26 和 27
pbmc <- subset(pbmc, subset = seurat_clusters %in% c("26", "27"), invert = TRUE)

# 2. 字典映射关系
cluster2celltype <- c(
  "0" = "Mature B cells", "1" = "NK cells", "2" = "Monocytes", 
  "3" = "Monocytes", "4" = "Naive CD4", "5" = "Effector CD8 cells", 
  "6" = "Monocytes", "7" = "Platelets", "8" = "Tregs", 
  "9" = "Naive CD8", "10" = "Mature B cells", "11" = "Effector CD8 cells", 
  "12" = "Monocytes", "13" = "Basophils", "14" = "Mature B cells", 
  "15" = "Exhausted CD8", "16" = "Naive Bcells", "17" = "Resting T cells", 
  "18" = "Tregs", "19" = "Erythrocytes", "20" = "Infiltrating Neutrophils", 
  "21" = "Dendritic cell", "22" = "Monocytes", "23" = "gamma delta", 
  "24" = "Activated Neutrophils", "25" = "gamma delta T17 cells", 
  "28" = "Epithelial cells"
)

# ==========================================
# 3. 合并相同 celltype 并重新顺序编号
# ==========================================

mapped_names <- unname(cluster2celltype[as.character(pbmc$seurat_clusters)])
pure_celltypes <- ifelse(is.na(mapped_names), as.character(pbmc$seurat_clusters), mapped_names)

unique_types <- sort(unique(pure_celltypes))
type2newid <- setNames(as.character(seq_along(unique_types)), unique_types)
id2type <- setNames(names(type2newid), type2newid)

pbmc$new_clusters <- unname(
  factor(
    type2newid[pure_celltypes],
    levels = as.character(1:length(unique_types))
  )
)

pbmc$celltype <- paste0(pbmc$new_clusters, ": ", pure_celltypes)

qsave(pbmc, 'pbmc_corrected.qs')
pbmc = qread('pbmc_corrected.qs')

# ==========================================
# 4. 进一步过滤和合并群
# ==========================================

# A. 剔除指定群
pbmc <- subset(pbmc, subset = new_clusters %in% c("1", "5", "6", "9", "13", "18"), invert = TRUE)

# B. 合并 4 和 15 为 CD8
pbmc$new_clusters <- as.character(pbmc$new_clusters)
pbmc$new_clusters[pbmc$new_clusters == "15"] <- "4"
id2type["4"] <- "CD8"

# ==========================================
# 5. 重新计算连续序号
# ==========================================

old_ids_present <- sort(as.numeric(unique(pbmc$new_clusters)))
old_ids_present <- as.character(old_ids_present)

current_types <- id2type[old_ids_present]

new_id_map <- setNames(as.character(seq_along(old_ids_present)), old_ids_present)

pbmc$new_clusters <- unname(
  factor(
    new_id_map[pbmc$new_clusters],
    levels = as.character(seq_along(old_ids_present))
  )
)

id2type <- setNames(current_types, as.character(seq_along(old_ids_present)))

pbmc$celltype <- paste0(pbmc$new_clusters, ": ", id2type[as.character(pbmc$new_clusters)])

pbmc$cell_annotation <- factor(
  unname(id2type[as.character(pbmc$new_clusters)]),
  levels = unname(id2type[as.character(seq_along(old_ids_present))])
)

qsave(pbmc, 'pbmc_recorrected.qs')

# ==========================================
# 准备替换图例文本的向量
# ==========================================

cluster_levels <- levels(pbmc$new_clusters)

legend_labels <- sapply(cluster_levels, function(x) {
  paste0(x, ": ", id2type[x])
})

# ==========================================
# 绘图: 基于合并后的 new_clusters 绘图
# ==========================================

p_total_anno <- DimPlot(pbmc, reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_discrete(labels = legend_labels) +
  ggtitle('Corrected PBMC') + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_blank()
  )

p_cold_anno <- DimPlot(
  subset(pbmc, subset = Group == "Cold_4C"),
  reduction = "umap",
  group.by = "new_clusters",
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("Cold_4C") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank()
  ) + 
  NoLegend()

p_rt_anno <- DimPlot(
  subset(pbmc, subset = Group == "RT_25C"),
  reduction = "umap",
  group.by = "new_clusters",
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("RT_25C") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank()
  ) + 
  NoLegend()

p_tn_anno <- DimPlot(
  subset(pbmc, subset = Group == "TN_30C"),
  reduction = "umap",
  group.by = "new_clusters",
  label = TRUE,
  repel = TRUE
) + 
  ggtitle("TN_30C") + 
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank()
  ) + 
  NoLegend()

p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

if (!dir.exists("pictures")) dir.create("pictures")

file_name_anno_png <- file.path("pictures", "UMAP_Grid_corrected_PBMC.png")
ggsave(filename = file_name_anno_png, plot = p_final_anno, width = 15, height = 11, dpi = 300)

file_name_anno_pdf <- file.path("pictures", "UMAP_Grid_corrected_PBMC_update1.pdf")
ggsave(filename = file_name_anno_pdf, plot = p_final_anno, width = 15, height = 11, dpi = 300)

# ==========================================
# 6. 保留原自动寻找 marker 并绘制 DotPlot
# ==========================================

pbmc = qread('pbmc_recorrected.qs')

Idents(pbmc) <- "cell_annotation"

message("Step 6: Finding all markers for each cluster automatically...")

pbmc.markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

top5_markers <- pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

top5_genes <- unique(top5_markers$gene)

write.csv(
  pbmc.markers,
  file = file.path("pictures", "PBMC_FindAllMarkers_all.csv"),
  row.names = FALSE
)

write.csv(
  top5_markers,
  file = file.path("pictures", "PBMC_FindAllMarkers_top5.csv"),
  row.names = FALSE
)

p_dotplot_auto <- DotPlot(pbmc, features = top5_genes, group.by = "cell_annotation") +
  scale_color_gradient(low = "white", high = "darkred") + 
  labs(
    size = "Fraction of cells\nin group (%)", 
    color = "Mean expression\nin group"
  ) +
  theme_bw() +
  ggtitle('PBMC DotPlot: Auto Top5 Markers') +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

file_name_dotplot_auto_png <- file.path("pictures", "DotPlot_PBMC_Top5_Markers_auto.png")
ggsave(filename = file_name_dotplot_auto_png, plot = p_dotplot_auto, width = 16, height = 6, dpi = 300)

file_name_dotplot_auto_pdf <- file.path("pictures", "DotPlot_PBMC_Top5_Markers_auto.pdf")
ggsave(filename = file_name_dotplot_auto_pdf, plot = p_dotplot_auto, width = 16, height = 6, dpi = 300)

message("Step 6 auto marker DotPlot saved.")

# ==========================================
# 7. 根据 pbmc_dotplot_correct.csv 映射表读取 marker 并绘制 DotPlot
#    DotPlot 的细胞类型名称使用 csv 表格中的名称
# ==========================================

message("Step 7: Reading marker mapping from pbmc_dotplot_correct.csv...")

marker_csv <- file.path(getwd(), "pbmc_dotplot_correct.csv")

if (!file.exists(marker_csv)) {
  stop("未找到 pbmc_dotplot_correct.csv，请把该文件放到当前工作目录。")
}

# ==========================================
# 7.1 读取 marker 映射表
# 要求 csv 至少包含两列：
# celltype：DotPlot 中希望显示的细胞类型名称
# markers：marker gene
# ==========================================

read_marker_csv <- function(path) {
  x <- tryCatch(
    {
      read.csv(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        fileEncoding = "UTF-8-BOM"
      )
    },
    error = function(e) {
      message("UTF-8-BOM 读取失败，尝试使用 GB18030 编码读取。")
      read.csv(
        path,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE,
        fileEncoding = "GB18030"
      )
    }
  )
  return(x)
}

marker_raw <- read_marker_csv(marker_csv)

colnames(marker_raw) <- trimws(colnames(marker_raw))

required_cols <- c("celltype", "markers")

if (!all(required_cols %in% colnames(marker_raw))) {
  stop(
    "pbmc_dotplot_correct.csv 必须包含列名：celltype 和 markers。\n",
    "当前检测到的列名为：",
    paste(colnames(marker_raw), collapse = ", ")
  )
}

marker_table <- marker_raw %>%
  select(
    marker_celltype = celltype,
    gene = markers
  ) %>%
  mutate(
    marker_celltype = trimws(as.character(marker_celltype)),
    gene = trimws(as.character(gene))
  ) %>%
  filter(
    !is.na(marker_celltype),
    !is.na(gene),
    marker_celltype != "",
    gene != ""
  )

# 按 csv 中出现顺序固定 DotPlot 细胞类型顺序
csv_celltype_order <- unique(marker_table$marker_celltype)

# 每个 csv 细胞类型只取前 5 个 marker
marker_table <- marker_table %>%
  mutate(marker_celltype = factor(marker_celltype, levels = csv_celltype_order)) %>%
  group_by(marker_celltype) %>%
  mutate(marker_rank = row_number()) %>%
  slice_head(n = 5) %>%
  ungroup() %>%
  arrange(marker_celltype, marker_rank)

# ==========================================
# 7.2 根据 csv 名称构建 DotPlot 专用细胞类型映射
# ==========================================
# 左侧：当前 pbmc$cell_annotation 中的名称
# 右侧：pbmc_dotplot_correct.csv 中希望显示的名称

annotation_to_table_name <- c(
  "Tregs" = "Treg",
  "Platelets" = "Platelet",
  "NK cells" = "NK",
  "Naive CD4" = "CD4",
  "Resting T cells" = "CD4",
  "Monocytes" = "Monocyte",
  "Mature B cells" = "B cell",
  "Naive Bcells" = "B cell",
  "Infiltrating Neutrophils" = "Neutrophil",
  "Activated Neutrophils" = "Neutrophil",
  "gamma delta" = "γδT cell",
  "gamma delta T17 cells" = "γδT cell",
  "Exhausted CD8" = "Exhausted CD8",
  "Dendritic cell" = "DC",
  "CD8" = "CD8",
  "Effector CD8 cells" = "CD8",
  "Naive CD8" = "CD8",
  "Basophils" = "Basophil"
)

# 如果 pbmc$cell_annotation 本身已经是 csv 里的名字，也允许直接匹配
exact_map <- setNames(csv_celltype_order, csv_celltype_order)
exact_map <- exact_map[setdiff(names(exact_map), names(annotation_to_table_name))]
annotation_to_table_name <- c(annotation_to_table_name, exact_map)

current_annotation <- as.character(pbmc$cell_annotation)

mapped_annotation <- unname(annotation_to_table_name[current_annotation])

unmapped_annotation <- unique(current_annotation[is.na(mapped_annotation)])

if (length(unmapped_annotation) > 0) {
  message("以下 pbmc$cell_annotation 没有设置到 csv 名称的映射：")
  print(unmapped_annotation)
  stop("请在 annotation_to_table_name 中补充这些细胞类型的映射关系。")
}

invalid_mapped_annotation <- setdiff(unique(mapped_annotation), csv_celltype_order)

if (length(invalid_mapped_annotation) > 0) {
  message("以下映射后的名称不在 pbmc_dotplot_correct.csv 的 celltype 列中：")
  print(invalid_mapped_annotation)
  stop("请检查 annotation_to_table_name 的右侧名称是否与 csv 中 celltype 完全一致。")
}

# 只保留当前 Seurat 对象中实际存在的 csv 细胞类型
dotplot_celltype_order <- csv_celltype_order[csv_celltype_order %in% unique(mapped_annotation)]

# 新增 DotPlot 专用分组列
# 这一步决定 DotPlot y 轴显示 csv 表格中的名称
pbmc$cell_annotation_dotplot <- factor(
  mapped_annotation,
  levels = dotplot_celltype_order
)

Idents(pbmc) <- "cell_annotation_dotplot"

message("DotPlot 使用的细胞类型顺序：")
print(levels(pbmc$cell_annotation_dotplot))

# ==========================================
# 7.3 只保留当前 DotPlot 中存在细胞类型对应的 marker
# ==========================================

marker_table_used <- marker_table %>%
  filter(marker_celltype %in% dotplot_celltype_order) %>%
  mutate(marker_celltype = factor(marker_celltype, levels = dotplot_celltype_order)) %>%
  arrange(marker_celltype, marker_rank)

marker_count_check <- marker_table_used %>%
  count(marker_celltype, name = "n_markers")

message("csv 中每个细胞类型用于绘图的 marker 数量：")
print(marker_count_check)

less_than_5 <- marker_count_check %>%
  filter(n_markers < 5)

if (nrow(less_than_5) > 0) {
  message("以下细胞类型少于 5 个 marker，请检查 csv：")
  print(less_than_5)
}

# ==========================================
# 7.4 基因名大小写修正
# ==========================================
# 例如 csv 里可能写 CD8a、CD14、themis，
# 但 Seurat 对象 rownames 里可能是 Cd8a、Cd14、Themis。
# 这里优先精确匹配，匹配不到再忽略大小写匹配。

match_gene_symbol <- function(input_genes, object_genes) {
  object_genes_lower <- tolower(object_genes)
  
  matched <- sapply(input_genes, function(g) {
    if (g %in% object_genes) {
      return(g)
    }
    
    idx <- match(tolower(g), object_genes_lower)
    
    if (is.na(idx)) {
      return(NA_character_)
    } else {
      return(object_genes[idx])
    }
  })
  
  unname(matched)
}

marker_table_used$gene_in_object <- match_gene_symbol(
  marker_table_used$gene,
  rownames(pbmc)
)

missing_genes <- marker_table_used %>%
  filter(is.na(gene_in_object)) %>%
  select(marker_celltype, gene)

if (nrow(missing_genes) > 0) {
  message("以下 marker 在 Seurat 对象 rownames(pbmc) 中未找到，已从 Step 7 DotPlot 中移除：")
  print(missing_genes)
}

marker_table_used <- marker_table_used %>%
  filter(!is.na(gene_in_object))

# 检查去除缺失基因后，每个细胞类型是否仍有 5 个 marker
marker_count_after_match <- marker_table_used %>%
  count(marker_celltype, name = "n_markers_after_match")

message("基因名匹配后每个细胞类型实际进入 DotPlot 的 marker 数量：")
print(marker_count_after_match)

less_than_5_after_match <- marker_count_after_match %>%
  filter(n_markers_after_match < 5)

if (nrow(less_than_5_after_match) > 0) {
  message("以下细胞类型在基因名匹配后少于 5 个 marker，通常是因为部分基因不在 rownames(pbmc) 中：")
  print(less_than_5_after_match)
}

# DotPlot 的 features 不能重复，同一个基因如果出现在多个细胞类型中，x 轴只会显示一次
duplicated_genes <- marker_table_used %>%
  count(gene_in_object, name = "n") %>%
  filter(n > 1)

if (nrow(duplicated_genes) > 0) {
  message("以下 marker 在多个细胞类型中重复出现，DotPlot 的 x 轴只会显示一次：")
  print(duplicated_genes)
}

top5_genes_csv <- unique(marker_table_used$gene_in_object)

if (length(top5_genes_csv) == 0) {
  stop("Step 7 没有可用于 DotPlot 的 marker genes，请检查 pbmc_dotplot_correct.csv 的细胞类型名称和基因名。")
}

message("Step 7 DotPlot 实际使用的 marker 数量：")
print(length(top5_genes_csv))

message("Step 7 DotPlot 使用的 marker genes:")
print(top5_genes_csv)

# 保存实际用于绘图的 marker 映射，方便核对
write.csv(
  marker_table_used,
  file = file.path("pictures", "DotPlot_marker_mapping_used_from_pbmc_dotplot_correct.csv"),
  row.names = FALSE
)

# ==========================================
# 7.5 绘制 DotPlot
# ==========================================

p_dotplot_csv <- DotPlot(
  pbmc,
  features = top5_genes_csv,
  group.by = "cell_annotation_dotplot"
) +
  scale_color_gradient(low = "white", high = "darkred") +
  labs(
    size = "Fraction of cells\nin group (%)",
    color = "Mean expression\nin group"
  ) +
  theme_bw() +
  ggtitle("PBMC DotPlot: Corrected Marker Genes") +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      face = "italic",
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.5),
    panel.grid.minor = element_blank()
  )

# ==========================================
# 7.6 保存图片
# ==========================================

plot_width <- max(16, length(top5_genes_csv) * 0.35)

file_name_dotplot_csv_png <- file.path("pictures", "DotPlot_PBMC_Top5_Markers_corrected_csv.png")
ggsave(
  filename = file_name_dotplot_csv_png,
  plot = p_dotplot_csv,
  width = plot_width,
  height = 6,
  dpi = 300
)

file_name_dotplot_csv_pdf <- file.path("pictures", "DotPlot_PBMC_Top5_Markers_corrected_csv.pdf")
ggsave(
  filename = file_name_dotplot_csv_pdf,
  plot = p_dotplot_csv,
  width = plot_width,
  height = 6,
  dpi = 300
)

message("Step 7 corrected csv marker DotPlot saved.")