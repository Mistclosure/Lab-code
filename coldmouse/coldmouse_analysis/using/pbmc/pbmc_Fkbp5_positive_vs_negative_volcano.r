# ==============================================================================
# PBMC Fkbp5+ vs Fkbp5- 火山图
#
# 比较对象：
#   Fkbp5+ cells vs Fkbp5- cells
#
# 输出：
#   1. 差异分析表
#   2. 火山图 png/pdf
#   3. 标注 top 10 差异基因
# ==============================================================================

library(Seurat)
library(qs)
library(dplyr)
library(ggplot2)
library(ggrepel)

# ==========================================
# 0. 参数区
# ==========================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')

input_qs <- "pbmc_recorrected.qs"
output_prefix <- "PBMC_Fkbp5_positive_vs_negative"

gene_target <- "Fkbp5"

# 判断 Fkbp5+ 的阈值
expr_threshold <- 0

# 可选：指定 assay
# 不指定则使用 Seurat 对象当前 DefaultAssay
assay_to_use <- NULL
# assay_to_use <- "RNA"

# 可选：只分析某些 group
# 不筛选则设为 NULL
group_col <- "Group"
groups_to_keep <- NULL
# groups_to_keep <- c("Cold_4C", "RT_25C", "TN_30C")

# 差异分析参数
min_pct <- 0.1
logfc_threshold_findmarkers <- 0

# 火山图显著性阈值
log2fc_cutoff <- 0.25
padj_cutoff <- 0.05

# 标注前 N 个基因
top_n_label <- 10

# 是否在 top10 标注中排除 Fkbp5 本身
# 因为 Fkbp5 是分组依据，本身通常一定显著，默认不标注它
exclude_target_gene_from_label <- TRUE

# 图片参数
plot_width <- 7
plot_height <- 6

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
# 3. 可选：按 Group 筛选
# ==========================================

if (!is.null(groups_to_keep)) {
  if (!group_col %in% colnames(obj@meta.data)) {
    stop("metadata 中没有分组列：", group_col)
  }
  
  obj <- subset(
    obj,
    subset = .data[[group_col]] %in% groups_to_keep
  )
  
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

# p_val_adj 如果出现 0，绘图时需要替换成极小值，避免 -log10(0) = Inf
min_nonzero_padj <- min(de_result$p_val_adj[de_result$p_val_adj > 0], na.rm = TRUE)

if (is.infinite(min_nonzero_padj)) {
  min_nonzero_padj <- .Machine$double.xmin
}

de_result <- de_result %>%
  mutate(
    log2FC = .data[[logfc_col]],
    p_val_adj_plot = ifelse(p_val_adj == 0, min_nonzero_padj * 0.1, p_val_adj),
    neg_log10_padj = -log10(p_val_adj_plot),
    change = case_when(
      p_val_adj < padj_cutoff & log2FC > log2fc_cutoff ~ paste0("Up in ", positive_label),
      p_val_adj < padj_cutoff & log2FC < -log2fc_cutoff ~ paste0("Up in ", negative_label),
      TRUE ~ "Not significant"
    )
  )

# 保存完整差异分析结果
write.csv(
  de_result,
  file = file.path("tables", paste0(output_prefix, "_FindMarkers.csv")),
  row.names = FALSE
)

# ==========================================
# 7. 选择前 10 个标注基因
# ==========================================

label_candidates <- de_result

if (exclude_target_gene_from_label) {
  label_candidates <- label_candidates %>%
    filter(tolower(gene) != tolower(gene_in_object))
}

label_candidates_sig <- label_candidates %>%
  filter(p_val_adj < padj_cutoff) %>%
  arrange(p_val_adj, desc(abs(log2FC)))

if (nrow(label_candidates_sig) >= top_n_label) {
  top_label_genes <- label_candidates_sig %>%
    slice_head(n = top_n_label)
} else {
  top_label_genes <- label_candidates %>%
    arrange(p_val_adj, desc(abs(log2FC))) %>%
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
# 8. 绘制火山图
# ==========================================

volcano_colors <- c(
  paste0("Up in ", positive_label) = "#D51F26",
  paste0("Up in ", negative_label) = "#3C5488FF",
  "Not significant" = "grey75"
)

p_volcano <- ggplot(
  de_result,
  aes(
    x = log2FC,
    y = neg_log10_padj,
    color = change
  )
) +
  geom_point(
    alpha = 0.75,
    size = 1.4
  ) +
  geom_vline(
    xintercept = c(-log2fc_cutoff, log2fc_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_text_repel(
    data = top_label_genes,
    aes(
      x = log2FC,
      y = neg_log10_padj,
      label = gene
    ),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    show.legend = FALSE
  ) +
  scale_color_manual(values = volcano_colors) +
  theme_classic(base_size = 13) +
  labs(
    title = paste0(positive_label, " vs ", negative_label),
    subtitle = paste0(
      "Top ", top_n_label, " genes labeled",
      ifelse(exclude_target_gene_from_label, paste0(" excluding ", gene_in_object), "")
    ),
    x = paste0("Average log2FC (", positive_label, " / ", negative_label, ")"),
    y = "-log10 adjusted p-value",
    color = NULL
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 17
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 11
    ),
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    axis.text = element_text(
      color = "black",
      size = 11
    ),
    legend.position = "top",
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
  )

# ==========================================
# 9. 保存图片
# ==========================================

file_name_png <- file.path(
  "pictures",
  paste0(output_prefix, "_volcano.png")
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
  paste0(output_prefix, "_volcano.pdf")
)

ggsave(
  filename = file_name_pdf,
  plot = p_volcano,
  width = plot_width,
  height = plot_height
)

message("火山图绘制完成：", file_name_png)