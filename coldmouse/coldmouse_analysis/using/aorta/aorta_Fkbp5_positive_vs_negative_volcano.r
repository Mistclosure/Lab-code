# ==============================================================================
# PBMC Fkbp5+ vs Fkbp5- 火山图
# 仿论文附图风格：
#   x 轴固定 -4 到 4
#   下方双向箭头从 0 出发
#   右侧：Enriched in Fkbp5+
#   左侧：Depleted in Fkbp5+
# ==============================================================================

library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(grid)

# ==========================================
# 0. 参数区
# ==========================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

input_qs <- "aorta_corrected.qs"
output_prefix <- "AORTA_Fkbp5_positive_vs_negative"

gene_target <- "Fkbp5"
expr_threshold <- 0

assay_to_use <- NULL
# assay_to_use <- "RNA"

group_col <- "Group"
groups_to_keep <- NULL
# groups_to_keep <- c("Cold_4C", "RT_25C", "TN_30C")

min_pct <- 0.1
logfc_threshold_findmarkers <- 0

log2fc_cutoff <- 0.25
padj_cutoff <- 0.05

top_n_label <- 10
exclude_target_gene_from_label <- TRUE

# 固定火山图横坐标范围
x_axis_limit <- 4

plot_width <- 8
plot_height <- 7

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

# ==========================================
# 3. 可选：筛选 Group
# ==========================================

if (!is.null(groups_to_keep)) {
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("metadata 中没有分组列：", group_col)
  }
  
  cells_keep <- rownames(obj@meta.data)[
    obj@meta.data[[group_col]] %in% groups_to_keep
  ]
  
  obj <- subset(obj, cells = cells_keep)
  
  message("已筛选指定 Group：")
  print(groups_to_keep)
  message("筛选后细胞数：", ncol(obj))
}

# ==========================================
# 4. 自动匹配 Fkbp5 基因名
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
# 5. 标记 Fkbp5+ / Fkbp5-
# ==========================================

expr_data <- FetchData(obj, vars = gene_in_object)

obj$target_gene_expr <- expr_data[[gene_in_object]]

positive_label <- paste0(gene_in_object, "+")
negative_label <- paste0(gene_in_object, "-")

obj$target_gene_status <- ifelse(
  obj$target_gene_expr > expr_threshold,
  positive_label,
  negative_label
)

obj$target_gene_status <- factor(
  obj$target_gene_status,
  levels = c(negative_label, positive_label)
)

message("Fkbp5+ / Fkbp5- 细胞数量：")
print(table(obj$target_gene_status))

if (sum(obj$target_gene_status == positive_label) == 0) {
  stop("没有筛选到 ", positive_label, " 细胞，请检查阈值或表达矩阵。")
}

if (sum(obj$target_gene_status == negative_label) == 0) {
  stop("没有筛选到 ", negative_label, " 细胞，请检查阈值或表达矩阵。")
}

# ==========================================
# 6. Fkbp5+ vs Fkbp5- 差异分析
# ==========================================

Idents(obj) <- "target_gene_status"

message("开始差异分析：", positive_label, " vs ", negative_label)

de_result <- FindMarkers(
  obj,
  ident.1 = positive_label,
  ident.2 = negative_label,
  only.pos = FALSE,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold_findmarkers
)

de_result$gene <- rownames(de_result)

# 兼容不同 Seurat 版本的 logFC 列名
if ("avg_log2FC" %in% colnames(de_result)) {
  logfc_col <- "avg_log2FC"
} else if ("avg_logFC" %in% colnames(de_result)) {
  logfc_col <- "avg_logFC"
} else {
  stop("FindMarkers 结果中未找到 avg_log2FC 或 avg_logFC 列。")
}

# 避免 p_val_adj = 0 导致 -log10(0) = Inf
min_nonzero_padj <- min(de_result$p_val_adj[de_result$p_val_adj > 0], na.rm = TRUE)

if (is.infinite(min_nonzero_padj)) {
  min_nonzero_padj <- .Machine$double.xmin
}

de_result <- de_result %>%
  mutate(
    log2FC = .data[[logfc_col]],
    p_val_adj_plot = ifelse(p_val_adj == 0, min_nonzero_padj * 0.1, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_plot),
    significant = p_val_adj < padj_cutoff & abs(log2FC) > log2fc_cutoff,
    significant_label = factor(
      ifelse(significant, "True", "False"),
      levels = c("True", "False")
    ),
    change = case_when(
      p_val_adj < padj_cutoff & log2FC > log2fc_cutoff ~ paste0("Enriched in ", positive_label),
      p_val_adj < padj_cutoff & log2FC < -log2fc_cutoff ~ paste0("Depleted in ", positive_label),
      TRUE ~ "Not significant"
    )
  )

write.csv(
  de_result,
  file = file.path("tables", paste0(output_prefix, "_FindMarkers.csv")),
  row.names = FALSE
)

# ==========================================
# 7. 选择火山图可见范围内的前 10 个标注基因
# ==========================================
# 注意：横坐标固定为 -4 到 4，所以标注也只从该范围内选择。
# 超出范围的基因仍保存在差异分析表里，但不显示在图中。

label_candidates <- de_result %>%
  filter(
    !is.na(log2FC),
    !is.na(neg_log10_padj),
    log2FC >= -x_axis_limit,
    log2FC <= x_axis_limit
  )

if (exclude_target_gene_from_label) {
  label_candidates <- label_candidates %>%
    filter(tolower(gene) != tolower(gene_in_object))
}

label_candidates_sig <- label_candidates %>%
  filter(significant) %>%
  arrange(p_val_adj, desc(abs(log2FC)))

if (nrow(label_candidates_sig) >= top_n_label) {
  top_label_genes <- label_candidates_sig %>%
    slice_head(n = top_n_label)
} else {
  top_label_genes <- bind_rows(
    label_candidates_sig,
    label_candidates %>%
      filter(!gene %in% label_candidates_sig$gene) %>%
      arrange(p_val_adj, desc(abs(log2FC)))
  ) %>%
    slice_head(n = top_n_label)
}

write.csv(
  top_label_genes,
  file = file.path("tables", paste0(output_prefix, "_top", top_n_label, "_label_genes.csv")),
  row.names = FALSE
)

message("火山图标注基因：")
print(top_label_genes[, c("gene", "log2FC", "p_val_adj", "change")])

# ==========================================
# 8. 准备绘图数据
# ==========================================

plot_data <- de_result %>%
  filter(
    !is.na(log2FC),
    !is.na(neg_log10_padj),
    log2FC >= -x_axis_limit,
    log2FC <= x_axis_limit
  )

message("由于横坐标固定为 -", x_axis_limit, " 到 ", x_axis_limit, "，以下数量的基因未显示在火山图中：")
message(nrow(de_result) - nrow(plot_data))

sig_colors <- setNames(
  c("red", "grey75"),
  c("True", "False")
)

# ==========================================
# 9. 火山图主体
# ==========================================

p_volcano_main <- ggplot(
  plot_data,
  aes(
    x = log2FC,
    y = neg_log10_padj
  )
) +
  geom_point(
    aes(color = significant_label),
    alpha = 0.8,
    size = 1.5
  ) +
  geom_text_repel(
    data = top_label_genes,
    aes(
      x = log2FC,
      y = neg_log10_padj,
      label = gene
    ),
    size = 4,
    fontface = "italic",
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.3,
    segment.color = "black",
    segment.linewidth = 0.4,
    min.segment.length = 0,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = sig_colors,
    breaks = c("True", "False"),
    labels = c("True", "False"),
    name = NULL
  ) +
  scale_x_continuous(
    limits = c(-x_axis_limit, x_axis_limit),
    breaks = seq(-x_axis_limit, x_axis_limit, by = 2),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = paste0(positive_label, " vs ", negative_label),
    x = NULL,
    y = expression(-log[10] * "(adjusted " * italic(P) * "-value)")
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.title.x = element_blank(),
    axis.text = element_text(
      color = "black",
      size = 12
    ),
    legend.position = c(0.10, 0.88),
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 12),
    plot.margin = margin(t = 20, r = 20, b = 5, l = 20)
  )

# ==========================================
# 10. 单独绘制底部箭头区域（修正版）
# ==========================================
# 目标：
# 1) 左右箭头从中心向两边发散
# 2) 中间留空放 log2FC
# 3) 左右文字下移，避免与箭头重叠

arrow_end <- x_axis_limit * 0.95   # 箭头终点
center_gap <- 0.75                 # 中间留白宽度，用于放 log2FC
arrow_y <- 0.78                    # 箭头所在高度
label_y <- 0.24                    # 左右说明文字高度
center_text_y <- 0.78              # log2FC 与箭头同高，但因为中间留白所以不会重叠

p_arrow <- ggplot() +
  # 左箭头：从中心左侧边界指向左端
  annotate(
    "segment",
    x = -center_gap,
    xend = -arrow_end,
    y = arrow_y,
    yend = arrow_y,
    arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  # 右箭头：从中心右侧边界指向右端
  annotate(
    "segment",
    x = center_gap,
    xend = arrow_end,
    y = arrow_y,
    yend = arrow_y,
    arrow = arrow(length = unit(0.20, "cm"), type = "closed"),
    linewidth = 0.6
  ) +
  # 中间的 log2FC（放在留白中）
  annotate(
    "text",
    x = 0,
    y = center_text_y,
    label = "log[2]*FC",
    parse = TRUE,
    size = 5.8,
    fontface = "bold"
  ) +
  # 左侧说明文字：下移，避免与箭头线重叠
  annotate(
    "text",
    x = -x_axis_limit + 0.12,
    y = label_y,
    label = paste0("Depleted in\n", positive_label),
    hjust = 0,
    vjust = 0.5,
    size = 4.9
  ) +
  # 右侧说明文字：下移，避免与箭头线重叠
  annotate(
    "text",
    x = x_axis_limit - 0.12,
    y = label_y,
    label = paste0("Enriched in\n", positive_label),
    hjust = 1,
    vjust = 0.5,
    size = 4.9
  ) +
  scale_x_continuous(
    limits = c(-x_axis_limit, x_axis_limit),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_void() +
  theme(
    plot.margin = margin(t = 0, r = 20, b = 20, l = 20)
  )
# ==========================================
# 11. 合并火山图主体和箭头区域
# ==========================================

p_volcano <- p_volcano_main / p_arrow +
  plot_layout(heights = c(1, 0.24))

# ==========================================
# 12. 保存图片
# ==========================================

file_name_png <- file.path(
  "pictures",
  paste0(output_prefix, "_volcano_paper_style.png")
)

ggsave(
  filename = file_name_png,
  plot = p_volcano,
  width = plot_width,
  height = plot_height,
  dpi = 300
)

file_name_pdf <- file.path(
  "pictures",
  paste0(output_prefix, "_volcano_paper_style.pdf")
)

ggsave(
  filename = file_name_pdf,
  plot = p_volcano,
  width = plot_width,
  height = plot_height
)

message("火山图绘制完成：", file_name_png)