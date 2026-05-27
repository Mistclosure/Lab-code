# ==============================================================================
# Script: Fkbp5+/- Monocyte 的 GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 通路打分与柱状图
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)
library(msigdbr)
library(ggplot2)

dir.create("pictures", showWarnings = FALSE)
dir.create("files", showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 读取保存好的 Monocytes qs 文件
# ------------------------------------------------------------------------------

print("🚀 步骤1: 读取保存好的单核细胞对象...")

mono_cells <- qread("pbmc_monocytes_sub-clustered.qs")

DefaultAssay(mono_cells) <- "RNA"

target_gene <- "Fkbp5"

if (!target_gene %in% rownames(mono_cells)) {
  stop("对象中找不到 Fkbp5，请检查基因名是否为 Fkbp5 或 FKBP5。")
}

print("  ✅ Monocytes 对象读取完成。")

# ------------------------------------------------------------------------------
# 2. 定义 Fkbp5+ / Fkbp5- Monocyte
# ------------------------------------------------------------------------------

print("🚀 步骤2: 根据 Fkbp5 counts 定义 Fkbp5+ / Fkbp5- Monocyte...")

# 兼容 Seurat v5 layer 和 Seurat v4 slot
get_counts_matrix <- function(obj, assay = "RNA") {
  tryCatch(
    GetAssayData(obj, assay = assay, layer = "counts"),
    error = function(e) {
      GetAssayData(obj, assay = assay, slot = "counts")
    }
  )
}

fkbp5_counts <- get_counts_matrix(mono_cells, assay = "RNA")[target_gene, ]

mono_cells$Fkbp5_Status <- ifelse(
  fkbp5_counts > 0,
  "Fkbp5+ Monocyte",
  "Fkbp5- Monocyte"
)

mono_cells$Fkbp5_Status <- factor(
  mono_cells$Fkbp5_Status,
  levels = c("Fkbp5- Monocyte", "Fkbp5+ Monocyte")
)

print(table(mono_cells$Fkbp5_Status))

print("  ✅ Fkbp5 状态定义完成。")

# ------------------------------------------------------------------------------
# 3. 读取 MSigDB mouse GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 通路基因
# ------------------------------------------------------------------------------

print("🚀 步骤3: 读取 GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 通路基因...")

gene_set_name <- "GOBP_ACTIVATION_OF_IMMUNE_RESPONSE"

# 兼容新版 / 旧版 msigdbr
msig_mouse_gobp <- tryCatch(
  {
    msigdbr(
      db_species = "MM",
      species = "Mus musculus",
      collection = "M5",
      subcollection = "GO:BP"
    )
  },
  error = function(e1) {
    message("新版 msigdbr 参数读取失败，尝试旧版 category/subcategory 参数...")
    msigdbr(
      species = "Mus musculus",
      category = "M5",
      subcategory = "GO:BP"
    )
  }
)

immune_genes <- msig_mouse_gobp %>%
  filter(gs_name == gene_set_name) %>%
  pull(gene_symbol) %>%
  unique()

if (length(immune_genes) == 0) {
  stop("没有从 msigdbr 中找到 GOBP_ACTIVATION_OF_IMMUNE_RESPONSE，请检查 msigdbr 版本或通路名。")
}

immune_genes_use <- intersect(immune_genes, rownames(mono_cells))

print(paste0("  通路总基因数: ", length(immune_genes)))
print(paste0("  在当前 Seurat 对象中匹配到的基因数: ", length(immune_genes_use)))

if (length(immune_genes_use) < 10) {
  stop("匹配到的通路基因少于 10 个，建议检查基因命名格式。")
}

write.csv(
  data.frame(gene = immune_genes_use),
  file = "files/GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_genes_used_in_monocytes.csv",
  row.names = FALSE
)

print("  ✅ 通路基因读取完成。")

# ------------------------------------------------------------------------------
# 4. 使用 AddModuleScore 计算免疫激活通路得分
# ------------------------------------------------------------------------------

print("🚀 步骤4: 使用 AddModuleScore 计算免疫通路得分...")

old_meta_cols <- colnames(mono_cells@meta.data)

set.seed(123)

mono_cells <- AddModuleScore(
  object = mono_cells,
  features = list(immune_genes_use),
  assay = "RNA",
  name = "GOBP_Activation_Immune_Response"
)

new_meta_cols <- setdiff(colnames(mono_cells@meta.data), old_meta_cols)
score_col <- new_meta_cols[grepl("GOBP_Activation_Immune_Response", new_meta_cols)]

if (length(score_col) != 1) {
  stop("无法唯一识别 AddModuleScore 生成的得分列，请检查 meta.data。")
}

print(paste0("  ✅ 得分列名: ", score_col))

# ------------------------------------------------------------------------------
# 5. 整理 Fkbp5+ / Fkbp5- Monocyte 的通路得分
# ------------------------------------------------------------------------------

print("🚀 步骤5: 整理绘图数据...")

plot_df <- mono_cells@meta.data %>%
  rownames_to_column("cell") %>%
  select(cell, Fkbp5_Status, all_of(score_col)) %>%
  rename(Immune_Activation_Score = all_of(score_col)) %>%
  filter(!is.na(Fkbp5_Status))

summary_df <- plot_df %>%
  group_by(Fkbp5_Status) %>%
  summarise(
    n = n(),
    mean_score = mean(Immune_Activation_Score, na.rm = TRUE),
    sd_score = sd(Immune_Activation_Score, na.rm = TRUE),
    sem_score = sd_score / sqrt(n),
    median_score = median(Immune_Activation_Score, na.rm = TRUE),
    .groups = "drop"
  )

write.csv(
  summary_df,
  file = "files/Summary_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.csv",
  row.names = FALSE
)

print(summary_df)

# Wilcoxon 检验
wilcox_res <- wilcox.test(
  Immune_Activation_Score ~ Fkbp5_Status,
  data = plot_df
)

p_label <- paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))

print(p_label)

# ------------------------------------------------------------------------------
# 6. 绘制箱线图
# 横坐标: Fkbp5+ / Fkbp5- Monocyte
# 纵坐标: GOBP_ACTIVATION_OF_IMMUNE_RESPONSE 通路得分
# ------------------------------------------------------------------------------

print("🚀 步骤6: 绘制免疫通路得分箱线图...")

# Wilcoxon 检验
wilcox_res <- wilcox.test(
  Immune_Activation_Score ~ Fkbp5_Status,
  data = plot_df
)

p_value <- wilcox_res$p.value

p_label <- case_when(
  p_value < 0.001 ~ "***",
  p_value < 0.01  ~ "**",
  p_value < 0.05  ~ "*",
  TRUE ~ "NS"
)

# 设置显著性标注位置
y_max <- max(plot_df$Immune_Activation_Score, na.rm = TRUE)
y_min <- min(plot_df$Immune_Activation_Score, na.rm = TRUE)
y_range <- y_max - y_min
y_text <- y_max + 0.12 * y_range
y_bracket <- y_max + 0.06 * y_range

# 红色虚线：这里用全部细胞的中位数作为参考线
# 如果你想改成 0，可以把 yintercept = median_score_ref 改成 yintercept = 0
median_score_ref <- median(plot_df$Immune_Activation_Score, na.rm = TRUE)

p_box <- ggplot(
  plot_df,
  aes(x = Fkbp5_Status, y = Immune_Activation_Score, fill = Fkbp5_Status)
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 16,
    outlier.size = 1.2,
    linewidth = 0.5,
    color = "black"
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.35,
    alpha = 0.35,
    color = "black"
  ) +
  geom_hline(
    yintercept = median_score_ref,
    linetype = "dashed",
    color = "red",
    linewidth = 0.5
  ) +
  # 显著性横线
  geom_segment(aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(aes(x = 1, xend = 1, y = y_bracket, yend = y_bracket - 0.03 * y_range),
               inherit.aes = FALSE, linewidth = 0.5) +
  geom_segment(aes(x = 2, xend = 2, y = y_bracket, yend = y_bracket - 0.03 * y_range),
               inherit.aes = FALSE, linewidth = 0.5) +
  annotate(
    "text",
    x = 1.5,
    y = y_text,
    label = p_label,
    size = 7,
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = c(
      "Fkbp5- Monocyte" = "grey80",
      "Fkbp5+ Monocyte" = "#E64B35"
    )
  ) +
  labs(
    title = "GOBP Activation of Immune Response",
    subtitle = paste0(
      "Fkbp5+ vs Fkbp5- Monocytes, Wilcoxon p = ",
      signif(p_value, 3)
    ),
    x = NULL,
    y = "Immune activation pathway score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1,
      face = "bold",
      size = 12
    ),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "none",
    plot.margin = margin(10, 15, 10, 10)
  )

ggsave(
  filename = "pictures/Boxplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.png",
  plot = p_box,
  width = 5.5,
  height = 5.5,
  dpi = 300
)

ggsave(
  filename = "pictures/Boxplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.pdf",
  plot = p_box,
  width = 5.5,
  height = 5.5
)

print("  ✅ 箱线图已保存。")
print("  pictures/Boxplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.png")
print("  pictures/Boxplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.pdf")
# ------------------------------------------------------------------------------
# 7. 保存带有通路得分的新 qs 对象
# ------------------------------------------------------------------------------

# qsave(
#   mono_cells,
#   "pbmc_monocytes_sub-clustered_Fkbp5_GOBP_ImmuneActivationScore.qs"
# )

print("✅ 全部完成！")
print("输出文件:")
print("  pictures/Barplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.png")
print("  pictures/Barplot_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.pdf")
print("  files/Summary_Fkbp5_Status_GOBP_ACTIVATION_OF_IMMUNE_RESPONSE_score.csv")
# print("  pbmc_monocytes_sub-clustered_Fkbp5_GOBP_ImmuneActivationScore.qs")