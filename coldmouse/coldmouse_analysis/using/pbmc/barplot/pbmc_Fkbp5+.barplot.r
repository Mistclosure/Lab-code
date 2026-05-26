# ==============================================================================
# 通用 Fkbp5+ 细胞组成百分比图脚本
# 横坐标：Percentage (%)
# 纵坐标：Group
# 填充：细胞类型注释
#
# 含义：
# 每个 Group 中，Fkbp5+ 细胞内部各细胞类型的组成百分比
# ==============================================================================

library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)

# ==========================================
# 0. 参数区：更换数据源时主要改这里
# ==========================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

input_qs <- "pbmc_recorrected.qs"
output_prefix <- "PBMC_Fkbp5_positive"

group_col <- "Group"
celltype_col <- "cell_annotation"

# 目标基因
gene_target <- "Fkbp5"

# 判断 Fkbp5+ 的阈值
# 通常归一化后的 data 层中 > 0 可认为表达
expr_threshold <- 0

# 可选：指定 assay，不指定则使用当前 DefaultAssay
assay_to_use <- NULL
# assay_to_use <- "RNA"

# 可选：手动指定 Group 顺序
#group_levels_manual <- NULL
group_levels_manual <- c('TN_30C','RT_25C','Cold_4C')

# 可选：手动指定细胞类型顺序
celltype_levels_manual <- NULL
# celltype_levels_manual <- c("Tregs", "Platelets", "NK cells", "Naive CD4", "Monocytes")

# 图片参数
plot_width <- 10
plot_height <- 5

# ==========================================
# 1. 创建输出目录
# ==========================================

if (!dir.exists("pictures")) dir.create("pictures")
if (!dir.exists("tables")) dir.create("tables")

# ==========================================
# 2. 读取 Seurat 对象
# ==========================================

if (!file.exists(input_qs)) {
  stop("未找到输入文件：", input_qs)
}

obj <- qread(input_qs)

if (!is.null(assay_to_use)) {
  if (!assay_to_use %in% Assays(obj)) {
    stop("对象中没有 assay：", assay_to_use)
  }
  DefaultAssay(obj) <- assay_to_use
}

message("使用数据源：", input_qs)
message("当前 DefaultAssay：", DefaultAssay(obj))
message("使用分组列：", group_col)
message("使用细胞注释列：", celltype_col)

# ==========================================
# 3. 检查 metadata 列
# ==========================================

if (!group_col %in% colnames(obj@meta.data)) {
  stop("metadata 中没有分组列：", group_col)
}

if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("metadata 中没有细胞注释列：", celltype_col)
}

# ==========================================
# 4. 匹配目标基因名
# ==========================================

match_gene_symbol <- function(input_gene, object_genes) {
  if (input_gene %in% object_genes) {
    return(input_gene)
  }
  
  idx <- match(tolower(input_gene), tolower(object_genes))
  
  if (is.na(idx)) {
    return(NA_character_)
  } else {
    return(object_genes[idx])
  }
}

gene_in_object <- match_gene_symbol(gene_target, rownames(obj))

if (is.na(gene_in_object)) {
  stop(
    "未在对象 rownames(obj) 中找到目标基因：", gene_target, "\n",
    "请检查基因名，例如 Fkbp5 / FKBP5 是否存在于表达矩阵。"
  )
}

message("目标基因匹配结果：", gene_target, " -> ", gene_in_object)

# ==========================================
# 5. 新增步骤：筛选 Fkbp5+ 细胞
# ==========================================

expr_data <- FetchData(obj, vars = gene_in_object)

obj$target_gene_expr <- expr_data[[gene_in_object]]

obj$target_gene_status <- ifelse(
  obj$target_gene_expr > expr_threshold,
  paste0(gene_in_object, "+"),
  paste0(gene_in_object, "-")
)

positive_label <- paste0(gene_in_object, "+")

positive_cells <- rownames(obj@meta.data)[obj$target_gene_status == positive_label]

message("筛选到 ", positive_label, " 细胞数量：", length(positive_cells))

if (length(positive_cells) == 0) {
  stop(
    "没有筛选到 ", positive_label, " 细胞。\n",
    "请检查 gene_target、expr_threshold 或表达矩阵。"
  )
}

# 后续所有计算都只基于 Fkbp5+ 细胞
obj_pos <- subset(obj, cells = positive_cells)

# ==========================================
# 6. 提取 Fkbp5+ 细胞 metadata
# ==========================================

meta_data <- obj_pos@meta.data %>%
  mutate(
    group = as.character(.data[[group_col]]),
    celltype = as.character(.data[[celltype_col]])
  ) %>%
  filter(
    !is.na(group),
    !is.na(celltype),
    group != "",
    celltype != ""
  )

# ==========================================
# 7. 设置 Group 顺序
# ==========================================

if (!is.null(group_levels_manual)) {
  group_levels <- group_levels_manual
  group_levels <- group_levels[group_levels %in% unique(meta_data$group)]
} else if (is.factor(obj@meta.data[[group_col]])) {
  group_levels <- levels(obj@meta.data[[group_col]])
  group_levels <- group_levels[group_levels %in% unique(meta_data$group)]
} else {
  group_levels <- unique(meta_data$group)
}

# ==========================================
# 8. 设置 celltype 顺序
# ==========================================

if (!is.null(celltype_levels_manual)) {
  celltype_levels <- celltype_levels_manual
  celltype_levels <- celltype_levels[celltype_levels %in% unique(meta_data$celltype)]
} else if (is.factor(obj@meta.data[[celltype_col]])) {
  celltype_levels <- levels(obj@meta.data[[celltype_col]])
  celltype_levels <- celltype_levels[celltype_levels %in% unique(meta_data$celltype)]
} else {
  celltype_levels <- unique(meta_data$celltype)
}

meta_data <- meta_data %>%
  mutate(
    group = factor(group, levels = group_levels),
    celltype = factor(celltype, levels = celltype_levels)
  )

# ==========================================
# 9. 计算每个 Group 内 Fkbp5+ 细胞的细胞类型比例
# ==========================================

cell_count <- meta_data %>%
  count(group, celltype, name = "cell_count")

group_total <- meta_data %>%
  count(group, name = "total_cells")

prop_data <- expand.grid(
  group = group_levels,
  celltype = celltype_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = factor(group, levels = group_levels),
    celltype = factor(celltype, levels = celltype_levels)
  ) %>%
  left_join(cell_count, by = c("group", "celltype")) %>%
  left_join(group_total, by = "group") %>%
  mutate(
    cell_count = ifelse(is.na(cell_count), 0, cell_count),
    total_cells = ifelse(is.na(total_cells), 0, total_cells),
    proportion = ifelse(total_cells > 0, cell_count / total_cells * 100, 0)
  )

write.csv(
  prop_data,
  file = file.path("tables", paste0(output_prefix, "_celltype_percentage_by_Group.csv")),
  row.names = FALSE
)

message("各 Group 中 ", positive_label, " 细胞数量：")
print(group_total)

message(positive_label, " 细胞内部各细胞类型比例：")
print(prop_data)

# ==========================================
# 10. 设置颜色
# ==========================================

celltype_colors <- scales::hue_pal()(length(celltype_levels))
names(celltype_colors) <- celltype_levels

# ==========================================
# 11. 绘制 Fkbp5+ 细胞水平堆叠百分比图
# ==========================================

p_positive <- ggplot(
  prop_data,
  aes(
    x = group,
    y = proportion,
    fill = celltype
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.65
  ) +
  coord_flip() +
  scale_fill_manual(values = celltype_colors) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  theme_classic() +
  labs(
    title = paste0("Cell Type Proportions in ", positive_label, " Cells"),
    x = NULL,
    y = paste0("Percentage of ", positive_label, " cells (%)"),
    fill = "Cell Type"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.text.y = element_text(
      size = 13,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 12,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "right",
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    legend.text = element_text(size = 10),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# ==========================================
# 12. 保存图片
# ==========================================

ggsave(
  filename = file.path("pictures", paste0(output_prefix, "_stacked.png")),
  plot = p_positive,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

ggsave(
  filename = file.path("pictures", paste0(output_prefix, "_stacked.pdf")),
  plot = p_positive,
  width = plot_width,
  height = plot_height
)

message(positive_label, " 细胞组成百分比图绘制完成。")