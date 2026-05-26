# ==========================================
# 07. PBMC 各细胞类型比例箱线图
# 横坐标：cell_annotation
# 纵坐标：Fraction of all cells (%)
# 分组：Group 中的所有分组
# ==========================================

library(Seurat)
library(qs)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

if (!dir.exists("pictures")) dir.create("pictures")
if (!dir.exists("tables")) dir.create("tables")

# ==========================================
# 1. 读取之前校正后的 PBMC 对象
# ==========================================

pbmc <- qread("pbmc_recorrected.qs")

# ==========================================
# 2. 设置列名
# ==========================================

group_col <- "Group"

# 你的前面代码里是 cell_annotation
# 如果你的对象里叫 cellannotation，也会自动兼容
if ("cell_annotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cell_annotation"
} else if ("cellannotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cellannotation"
} else {
  stop("没有找到 cell_annotation 或 cellannotation 列，请检查 pbmc@meta.data。")
}

if (!group_col %in% colnames(pbmc@meta.data)) {
  stop("没有找到 Group 列，请检查 pbmc@meta.data。")
}

# ==========================================
# 3. 自动识别样本列
# ==========================================
# 箱线图需要每个 Group 下有多个样本。
# 这里优先寻找常见样本列。
# 如果你的样本列不是这些名字，请手动修改 sample_col。

candidate_sample_cols <- c(
  "sample", "Sample", "sample_id", "SampleID",
  "orig.ident", "orig_ident",
  "mouse", "Mouse",
  "donor", "Donor",
  "batch", "Batch"
)

sample_col <- candidate_sample_cols[
  candidate_sample_cols %in% colnames(pbmc@meta.data)
][1]

if (is.na(sample_col)) {
  message("没有自动识别到样本列，将使用 Group 作为 sample_id。")
  message("注意：如果每个 Group 只有一个值，箱线图不会有真实的组内分布。")
  sample_col <- group_col
}

message("用于计算比例的细胞注释列：", celltype_col)
message("用于分组展示的列：", group_col)
message("用于样本重复的列：", sample_col)

# ==========================================
# 4. 提取 metadata
# ==========================================

meta_df <- pbmc@meta.data %>%
  mutate(
    cell_id = rownames(pbmc@meta.data),
    celltype = as.character(.data[[celltype_col]]),
    group = as.character(.data[[group_col]]),
    sample_id = as.character(.data[[sample_col]])
  ) %>%
  filter(
    !is.na(celltype),
    !is.na(group),
    !is.na(sample_id),
    celltype != "",
    group != "",
    sample_id != ""
  )

# 保持 cell_annotation 原来的顺序
if (is.factor(pbmc@meta.data[[celltype_col]])) {
  celltype_order <- levels(pbmc@meta.data[[celltype_col]])
  celltype_order <- celltype_order[celltype_order %in% unique(meta_df$celltype)]
} else {
  celltype_order <- unique(meta_df$celltype)
}

# 保持 Group 原来的顺序
if (is.factor(pbmc@meta.data[[group_col]])) {
  group_order <- levels(pbmc@meta.data[[group_col]])
  group_order <- group_order[group_order %in% unique(meta_df$group)]
} else {
  group_order <- unique(meta_df$group)
}

meta_df <- meta_df %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    group = factor(group, levels = group_order)
  )

# ==========================================
# 5. 计算每个样本中各细胞类型占全部 PBMC 细胞的比例
# ==========================================

sample_info <- meta_df %>%
  distinct(sample_id, group)

sample_total <- meta_df %>%
  count(sample_id, name = "n_total")

celltype_count <- meta_df %>%
  count(sample_id, group, celltype, name = "n_cell")

# 补全 0 值：
# 如果某个样本没有某类细胞，也保留为 0
plot_df <- expand.grid(
  sample_id = unique(meta_df$sample_id),
  celltype = celltype_order,
  stringsAsFactors = FALSE
) %>%
  left_join(sample_info, by = "sample_id") %>%
  left_join(celltype_count, by = c("sample_id", "group", "celltype")) %>%
  left_join(sample_total, by = "sample_id") %>%
  mutate(
    n_cell = ifelse(is.na(n_cell), 0, n_cell),
    fraction = n_cell / n_total,
    fraction_percent = fraction * 100,
    celltype = factor(celltype, levels = celltype_order),
    group = factor(group, levels = group_order)
  )

# 保存比例表，方便后续统计或检查
write.csv(
  plot_df,
  file = file.path("tables", "PBMC_cell_fraction_by_sample.csv"),
  row.names = FALSE
)

message("每个 Group 的样本数量：")
print(
  plot_df %>%
    distinct(sample_id, group) %>%
    count(group, name = "n_samples")
)

# ==========================================
# 6. 绘制箱线图
# ==========================================

p_fraction_boxplot <- ggplot(
  plot_df,
  aes(
    x = celltype,
    y = fraction_percent,
    fill = group
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_point(
    aes(color = group),
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width = 0.8
    ),
    size = 1.8,
    alpha = 0.85
  ) +
  labs(
    title = "PBMC Cell Fraction",
    x = NULL,
    y = "Fraction of all cells (%)",
    fill = "Group",
    color = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 18,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ==========================================
# 7. 保存图片
# ==========================================

plot_width <- max(10, length(celltype_order) * 0.8)

ggsave(
  filename = file.path("pictures", "PBMC_cell_fraction_boxplot_by_Group.png"),
  plot = p_fraction_boxplot,
  width = plot_width,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path("pictures", "PBMC_cell_fraction_boxplot_by_Group.pdf"),
  plot = p_fraction_boxplot,
  width = plot_width,
  height = 7
)

message("PBMC cell fraction boxplot saved.")






# ==========================================
# PBMC 各细胞类型比例图
# 横坐标：cell_annotation
# 纵坐标：Fraction of all cells (%)
# 展示：Group 中所有分组
# 适用于每个 Group 只有一个整体样本的情况
# ==========================================

library(Seurat)
library(qs)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

if (!dir.exists("pictures")) dir.create("pictures")
if (!dir.exists("tables")) dir.create("tables")

# ==========================================
# 1. 读取 PBMC 对象
# ==========================================

pbmc <- qread("pbmc_recorrected.qs")

# ==========================================
# 2. 检查必要列
# ==========================================

if (!"Group" %in% colnames(pbmc@meta.data)) {
  stop("pbmc@meta.data 中没有 Group 列。")
}

if ("cell_annotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cell_annotation"
} else if ("cellannotation" %in% colnames(pbmc@meta.data)) {
  celltype_col <- "cellannotation"
} else {
  stop("pbmc@meta.data 中没有 cell_annotation 或 cellannotation 列。")
}

# ==========================================
# 3. 提取 metadata
# ==========================================

meta_df <- pbmc@meta.data %>%
  mutate(
    celltype = as.character(.data[[celltype_col]]),
    group = as.character(Group)
  ) %>%
  filter(
    !is.na(celltype),
    !is.na(group),
    celltype != "",
    group != ""
  )

# 保持细胞类型原有顺序
if (is.factor(pbmc@meta.data[[celltype_col]])) {
  celltype_order <- levels(pbmc@meta.data[[celltype_col]])
  celltype_order <- celltype_order[celltype_order %in% unique(meta_df$celltype)]
} else {
  celltype_order <- unique(meta_df$celltype)
}

# 保持 Group 原有顺序
if (is.factor(pbmc@meta.data$Group)) {
  group_order <- levels(pbmc@meta.data$Group)
  group_order <- group_order[group_order %in% unique(meta_df$group)]
} else {
  group_order <- unique(meta_df$group)
}

meta_df <- meta_df %>%
  mutate(
    celltype = factor(celltype, levels = celltype_order),
    group = factor(group, levels = group_order)
  )

# ==========================================
# 4. 按 Group 计算各细胞类型占比
# ==========================================

group_total <- meta_df %>%
  count(group, name = "n_total")

cell_count <- meta_df %>%
  count(group, celltype, name = "n_cell")

plot_df <- expand.grid(
  group = group_order,
  celltype = celltype_order,
  stringsAsFactors = FALSE
) %>%
  mutate(
    group = factor(group, levels = group_order),
    celltype = factor(celltype, levels = celltype_order)
  ) %>%
  left_join(cell_count, by = c("group", "celltype")) %>%
  left_join(group_total, by = "group") %>%
  mutate(
    n_cell = ifelse(is.na(n_cell), 0, n_cell),
    fraction = n_cell / n_total,
    fraction_percent = fraction * 100
  )

write.csv(
  plot_df,
  file = file.path("tables", "PBMC_cell_fraction_by_Group.csv"),
  row.names = FALSE
)

# ==========================================
# 5. 绘制分组柱状图
# ==========================================

p_fraction_bar <- ggplot(
  plot_df,
  aes(
    x = celltype,
    y = fraction_percent,
    fill = group
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.2
  ) +
  labs(
    title = "PBMC Cell Fraction",
    x = NULL,
    y = "Fraction of all cells (%)",
    fill = "Group"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      size = 20,
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

plot_width <- max(10, length(celltype_order) * 0.8)

ggsave(
  filename = file.path("pictures", "PBMC_cell_fraction_barplot_by_Group.png"),
  plot = p_fraction_bar,
  width = plot_width,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path("pictures", "PBMC_cell_fraction_barplot_by_Group.pdf"),
  plot = p_fraction_bar,
  width = plot_width,
  height = 7
)

message("PBMC cell fraction barplot saved.")