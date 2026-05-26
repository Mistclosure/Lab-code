# ==============================================================================
# 不同细胞类型中 Fkbp5+ / Fkbp5- 细胞比例图
#
# 含义：
# 每个 Group 中，每一种 cell_annotation 内部，
# Fkbp5+ 和 Fkbp5- 细胞分别占该细胞类型的百分比。
#
# 图形：
# 横坐标：Percentage (%)
# 纵坐标：cell_annotation
# 分面：Group，例如 Cold_4C / RT_25C / TN_30C
# 图例：Fkbp5+ / Fkbp5-
# ==============================================================================

library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)

# ==========================================
# 0. 参数区
# ==========================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

input_qs <- "aorta_corrected.qs"
output_prefix <- "AORTA"

group_col <- "Group"
celltype_col <- "cell_annotation"

gene_target <- "Fkbp5"
expr_threshold <- 0

# 如果不指定 assay，则使用对象当前 DefaultAssay
assay_to_use <- NULL
# assay_to_use <- "RNA"

# 手动指定 Group 顺序
group_levels_manual <- c("Cold_4C", "RT_25C", "TN_30C")

# 如果想手动指定细胞类型顺序，可以在这里写
# 不想指定则设为 NULL，默认使用 cell_annotation 原来的 factor 顺序
celltype_levels_manual <- NULL

# 图片大小
plot_width <- 12
plot_height <- 12

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
# 4. 自动匹配目标基因名
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
    "请检查基因名是否为 Fkbp5 / FKBP5，或检查 DefaultAssay 是否正确。"
  )
}

message("目标基因匹配结果：", gene_target, " -> ", gene_in_object)

# ==========================================
# 5. 提取 Fkbp5 表达并标记 Fkbp5+ / Fkbp5-
# ==========================================

expr_data <- FetchData(obj, vars = gene_in_object)

obj$target_gene_expr <- expr_data[[gene_in_object]]

obj$target_gene_status <- ifelse(
  obj$target_gene_expr > expr_threshold,
  paste0(gene_in_object, "+"),
  paste0(gene_in_object, "-")
)

obj$target_gene_status <- factor(
  obj$target_gene_status,
  levels = c(paste0(gene_in_object, "-"), paste0(gene_in_object, "+"))
)

# ==========================================
# 6. 提取 metadata
# ==========================================

meta_data <- obj@meta.data %>%
  mutate(
    group = as.character(.data[[group_col]]),
    celltype = as.character(.data[[celltype_col]]),
    fkbp5_status = as.character(target_gene_status)
  ) %>%
  filter(
    !is.na(group),
    !is.na(celltype),
    !is.na(fkbp5_status),
    group != "",
    celltype != "",
    fkbp5_status != ""
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

status_levels <- c(paste0(gene_in_object, "-"), paste0(gene_in_object, "+"))

meta_data <- meta_data %>%
  mutate(
    group = factor(group, levels = group_levels),
    celltype = factor(celltype, levels = celltype_levels),
    fkbp5_status = factor(fkbp5_status, levels = status_levels)
  )

# ==========================================
# 9. 计算每个 Group + celltype 内 Fkbp5+ / Fkbp5- 比例
# ==========================================

count_data <- meta_data %>%
  count(group, celltype, fkbp5_status, name = "cell_count")

total_data <- meta_data %>%
  count(group, celltype, name = "total_cells")

plot_data <- expand.grid(
  group = group_levels,
  celltype = celltype_levels,
  fkbp5_status = status_levels,
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = factor(group, levels = group_levels),
    celltype = factor(celltype, levels = celltype_levels),
    fkbp5_status = factor(fkbp5_status, levels = status_levels)
  ) %>%
  left_join(count_data, by = c("group", "celltype", "fkbp5_status")) %>%
  left_join(total_data, by = c("group", "celltype")) %>%
  mutate(
    cell_count = ifelse(is.na(cell_count), 0, cell_count),
    total_cells = ifelse(is.na(total_cells), 0, total_cells),
    percentage = ifelse(total_cells > 0, cell_count / total_cells * 100, 0)
  )

# 保存统计表
write.csv(
  plot_data,
  file = file.path(
    "tables",
    paste0(output_prefix, "_", gene_in_object, "_positive_negative_percentage_by_celltype_and_Group.csv")
  ),
  row.names = FALSE
)

message("每个 Group + celltype 中 Fkbp5+ / Fkbp5- 比例：")
print(plot_data)

# 单独输出 Fkbp5+ 比例表，方便后续做统计或汇报
positive_ratio_table <- plot_data %>%
  filter(fkbp5_status == paste0(gene_in_object, "+")) %>%
  select(group, celltype, positive_cells = cell_count, total_cells, positive_percentage = percentage)

write.csv(
  positive_ratio_table,
  file = file.path(
    "tables",
    paste0(output_prefix, "_", gene_in_object, "_positive_ratio_by_celltype_and_Group.csv")
  ),
  row.names = FALSE
)

message("Fkbp5+ 比例表：")
print(positive_ratio_table)

# ==========================================
# 10. 设置颜色
# ==========================================

status_levels <- c(
  paste0(gene_in_object, "-"),
  paste0(gene_in_object, "+")
)

status_colors <- setNames(
  c("grey80", "#D51F26"),
  status_levels
)

# ==========================================
# 11. 绘制水平堆叠百分比图
# ==========================================

p_fkbp5_ratio <- ggplot(
  plot_data,
  aes(
    x = celltype,
    y = percentage,
    fill = fkbp5_status
  )
) +
  geom_col(
    color = "black",
    linewidth = 0.25,
    width = 0.7
  ) +
  coord_flip() +
  facet_wrap(~ group, ncol = 1) +
  scale_fill_manual(
    values = status_colors,
    breaks = c(paste0(gene_in_object, "+"), paste0(gene_in_object, "-")),
    labels = c(paste0(gene_in_object, "+"), paste0(gene_in_object, "-"))
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, by = 20),
    expand = c(0, 0)
  ) +
  theme_classic() +
  labs(
    title = paste0(gene_in_object, "+ Cell Proportion in Each Cell Type"),
    x = NULL,
    y = "Percentage within each cell type (%)",
    fill = NULL
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 17
    ),
    strip.text = element_text(
      size = 13,
      face = "bold"
    ),
    strip.background = element_rect(
      fill = "grey90",
      color = "black",
      linewidth = 0.3
    ),
    axis.text.y = element_text(
      size = 11,
      face = "bold",
      color = "black"
    ),
    axis.text.x = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "top",
    legend.text = element_text(
      size = 12,
      face = "bold"
    ),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# ==========================================
# 12. 保存图片
# ==========================================

file_name_png <- file.path(
  "pictures",
  paste0(output_prefix, "_", gene_in_object, "_positive_negative_percentage_by_celltype_and_Group.png")
)

ggsave(
  filename = file_name_png,
  plot = p_fkbp5_ratio,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

file_name_pdf <- file.path(
  "pictures",
  paste0(output_prefix, "_", gene_in_object, "_positive_negative_percentage_by_celltype_and_Group.pdf")
)

ggsave(
  filename = file_name_pdf,
  plot = p_fkbp5_ratio,
  width = plot_width,
  height = plot_height
)

message("Fkbp5+ / Fkbp5- 比例图绘制完成：", file_name_png)