# ==============================================================================
# 通用细胞组成百分比图脚本
# 横坐标：Percentage (%)
# 纵坐标：Group
# 填充：细胞类型注释
# ==============================================================================

library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)

# ==========================================
# 0. 参数区：更换数据源时主要改这里
# ==========================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

input_qs <- "aorta_corrected.qs"
output_prefix <- "AORTA"

group_col <- "Group"
celltype_col <- "cell_annotation"

# 可选：手动指定 group 顺序
# 如果不想指定，设为 NULL
group_levels_manual <- NULL
# group_levels_manual <- c("RT_25C", "Cold_4C", "TN_30C")

# 可选：手动指定细胞类型顺序
# 如果不想指定，设为 NULL
celltype_levels_manual <- NULL
# celltype_levels_manual <- c("Tregs", "Platelets", "NK cells", "Naive CD4", "Monocytes")

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
obj$Group = factor(obj$Group,levels = c('TN_30C','RT_25C','Cold_4C'))
# ==========================================
# 3. 检查 metadata 列
# ==========================================

if (!group_col %in% colnames(obj@meta.data)) {
  stop("metadata 中没有分组列：", group_col)
}

if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("metadata 中没有细胞注释列：", celltype_col)
}

message("使用数据源：", input_qs)
message("使用分组列：", group_col)
message("使用细胞注释列：", celltype_col)

# ==========================================
# 4. 提取 metadata
# ==========================================

meta_data <- obj@meta.data %>%
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
# 5. 设置 Group 顺序
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
# 6. 设置 celltype 顺序
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
# 7. 计算每个 Group 内各细胞类型比例
# ==========================================

prop_data <- meta_data %>%
  group_by(group, celltype) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells * 100
  ) %>%
  ungroup()

# 补全缺失组合
prop_data <- prop_data %>%
  complete(
    group = factor(group_levels, levels = group_levels),
    celltype = factor(celltype_levels, levels = celltype_levels),
    fill = list(
      cell_count = 0,
      total_cells = NA,
      proportion = 0
    )
  ) %>%
  group_by(group) %>%
  mutate(
    total_cells = ifelse(
      is.na(total_cells),
      sum(cell_count),
      total_cells
    )
  ) %>%
  ungroup()

# 保存比例表
write.csv(
  prop_data,
  file = file.path("tables", paste0(output_prefix, "_celltype_percentage_by_Group.csv")),
  row.names = FALSE
)

message("各 Group 的细胞总数：")
print(
  meta_data %>%
    count(group, name = "total_cells")
)

# ==========================================
# 8. 设置颜色
# ==========================================

celltype_colors <- scales::hue_pal()(length(celltype_levels))
names(celltype_colors) <- celltype_levels

# ==========================================
# 9. 绘制水平堆叠百分比图
# ==========================================

p_percentage <- ggplot(
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
    title = paste0(output_prefix, " Cell Type Proportions"),
    x = NULL,
    y = "Percentage of cells (%)",
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
# 10. 保存图片
# ==========================================

ggsave(
  filename = file.path("pictures", paste0(output_prefix, "_celltype_percentage_by_Group_stacked.png")),
  plot = p_percentage,
  width = 10,
  height = 5,
  dpi = 300
)

ggsave(
  filename = file.path("pictures", paste0(output_prefix, "_celltype_percentage_by_Group_stacked.pdf")),
  plot = p_percentage,
  width = 10,
  height = 5
)

message("比例图绘制完成。")