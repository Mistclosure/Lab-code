# ==============================================================================
# PBMC 细胞类型组成比例图
# 横坐标：Percentage (%)
# 纵坐标：Group
# 填充：cell_annotation
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

library(Seurat)
library(qs)
library(tidyverse)
library(ggplot2)

print("🚀 开始计算并绘制 PBMC 各 Group 的细胞组成百分比图...")

# ==========================================
# 1. 读取前面校正后的 PBMC 对象
# ==========================================

if (file.exists("pbmc_recorrected.qs")) {
  pbmc <- qread("pbmc_recorrected.qs")
} else {
  stop("❌ 未找到 pbmc_recorrected.qs，请先运行前面的 PBMC 注释校正脚本。")
}
pbmc$Group = factor(pbmc$Group,levels = c('TN_30C','RT_25C','Cold_4C'))
# ==========================================
# 2. 检查 metadata 列
# ==========================================

if (!"Group" %in% colnames(pbmc@meta.data)) {
  stop("❌ pbmc@meta.data 中没有 Group 列，请检查 metadata。")
}

if ("cell_annotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cell_annotation"
} else if ("cellannotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cellannotation"
} else {
  stop("❌ pbmc@meta.data 中没有 cell_annotation 或 cellannotation 列。")
}

print(paste("✅ 使用细胞注释列:", celltype_col))
print("✅ 使用分组列: Group")

# ==========================================
# 3. 提取 metadata 并整理
# ==========================================

meta_data <- pbmc@meta.data %>%
  mutate(
    Group = as.character(Group),
    cell_annotation_plot = as.character(.data[[celltype_col]])
  ) %>%
  filter(
    !is.na(Group),
    !is.na(cell_annotation_plot),
    Group != "",
    cell_annotation_plot != ""
  )

# 保持 Group 原有顺序
if (is.factor(pbmc@meta.data$Group)) {
  group_levels <- levels(pbmc@meta.data$Group)
  group_levels <- group_levels[group_levels %in% unique(meta_data$Group)]
} else {
  group_levels <- unique(meta_data$Group)
}

# 保持 cell_annotation 原有顺序
if (is.factor(pbmc@meta.data[[celltype_col]])) {
  celltype_levels <- levels(pbmc@meta.data[[celltype_col]])
  celltype_levels <- celltype_levels[celltype_levels %in% unique(meta_data$cell_annotation_plot)]
} else {
  celltype_levels <- unique(meta_data$cell_annotation_plot)
}

meta_data <- meta_data %>%
  mutate(
    Group = factor(Group, levels = group_levels),
    cell_annotation_plot = factor(cell_annotation_plot, levels = celltype_levels)
  )

# ==========================================
# 4. 计算每个 Group 内各细胞类型的百分比
# ==========================================

prop_data <- meta_data %>%
  group_by(Group, cell_annotation_plot) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(
    total_cells = sum(cell_count),
    proportion = cell_count / total_cells * 100
  ) %>%
  ungroup()

# 补全缺失组合：
# 如果某个 Group 中没有某类细胞，比例记为 0
prop_data <- prop_data %>%
  complete(
    Group = factor(group_levels, levels = group_levels),
    cell_annotation_plot = factor(celltype_levels, levels = celltype_levels),
    fill = list(
      cell_count = 0,
      total_cells = NA,
      proportion = 0
    )
  ) %>%
  group_by(Group) %>%
  mutate(
    total_cells = ifelse(
      is.na(total_cells),
      sum(cell_count),
      total_cells
    )
  ) %>%
  ungroup()

# 保存比例表
if (!dir.exists("tables")) dir.create("tables")

write.csv(
  prop_data,
  file = file.path("tables", "PBMC_cell_annotation_percentage_by_Group.csv"),
  row.names = FALSE
)

print("📊 各 Group 的细胞总数：")
print(
  meta_data %>%
    count(Group, name = "total_cells")
)

print("📊 各 Group 内各细胞类型比例：")
print(prop_data)

# ==========================================
# 5. 设置颜色
# ==========================================

celltype_colors <- scales::hue_pal()(length(celltype_levels))
names(celltype_colors) <- celltype_levels

# ==========================================
# 6. 绘制水平堆叠百分比图
# ==========================================

p_pbmc_percentage <- ggplot(
  prop_data,
  aes(
    x = Group,
    y = proportion,
    fill = cell_annotation_plot
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
    title = "PBMC Cell Type Proportions",
    x = NULL,
    y = "Percentage of PBMC cells (%)",
    fill = "Cell Annotation"
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
# 7. 保存图片
# ==========================================

if (!dir.exists("pictures")) dir.create("pictures")

file_name_png <- file.path("pictures", "PBMC_cell_annotation_percentage_by_Group_stacked.png")
ggsave(
  filename = file_name_png,
  plot = p_pbmc_percentage,
  width = 10,
  height = 5,
  dpi = 300
)

file_name_pdf <- file.path("pictures", "PBMC_cell_annotation_percentage_by_Group_stacked.pdf")
ggsave(
  filename = file_name_pdf,
  plot = p_pbmc_percentage,
  width = 10,
  height = 5
)

print(paste("✅ PBMC 细胞组成百分比图已保存至:", file_name_png))