# ==============================================================================
# Script: Aorta new_clusters Cold vs RT DEG
# 比较:
#   Cold cluster 7  vs RT cluster 7
#   Cold cluster 13 vs RT cluster 13
#   Cold cluster 16 vs RT cluster 16
#   Cold cluster 7+13 vs RT cluster 7+13
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)

if (!dir.exists("files")) dir.create("files")

# ------------------------------------------------------------------------------
# 0. 参数设置
# ------------------------------------------------------------------------------

group_cold <- "Cold_4C"
group_rt   <- "RT_25C"

fc_threshold <- 0.25
p_val_threshold <- 0.05

# ------------------------------------------------------------------------------
# 1. 读取 Aorta 对象
# ------------------------------------------------------------------------------

print("🚀 读取 Aorta 对象...")

aorta_obj <- qread("aorta_corrected.qs")
DefaultAssay(aorta_obj) <- "RNA"

if (!"Group" %in% colnames(aorta_obj@meta.data)) {
  stop("metadata 中没有 Group 列，请检查对象。")
}

if (!"new_clusters" %in% colnames(aorta_obj@meta.data)) {
  stop("metadata 中没有 new_clusters 列，请检查对象。")
}

# 只保留 Cold 和 RT
aorta_obj <- subset(aorta_obj, subset = Group %in% c(group_cold, group_rt))

# 确保 new_clusters 是字符型，避免 7 和 "7" 匹配问题
aorta_obj$new_clusters <- as.character(aorta_obj$new_clusters)

print("✅ 对象读取完成")
print(table(aorta_obj$Group))
print(table(aorta_obj$new_clusters, aorta_obj$Group))

# ------------------------------------------------------------------------------
# 2. 定义 DEG 分析函数
# ------------------------------------------------------------------------------

run_cluster_deg <- function(
    obj,
    clusters_use,
    comparison_name,
    group_cold = "Cold_4C",
    group_rt = "RT_25C",
    fc_threshold = 0.25,
    p_val_threshold = 0.05
) {
  
  print(paste0("🚀 开始分析: ", comparison_name))
  
  # 提取指定 cluster
  sub_obj <- subset(
    obj,
    subset = new_clusters %in% clusters_use
  )
  
  # 检查分组是否存在
  group_table <- table(sub_obj$Group)
  print(group_table)
  
  if (!all(c(group_cold, group_rt) %in% names(group_table))) {
    stop(paste0(
      comparison_name,
      " 中缺少 Cold 或 RT 分组，请检查该 cluster 是否同时存在两个组。"
    ))
  }
  
  Idents(sub_obj) <- "Group"
  
  # 差异分析：Cold vs RT
  deg_results <- FindMarkers(
    sub_obj,
    ident.1 = group_cold,
    ident.2 = group_rt,
    logfc.threshold = 0,
    min.pct = 0.1,
    verbose = FALSE
  )
  
  deg_results$gene <- rownames(deg_results)
  
  deg_results <- deg_results %>%
    mutate(
      Significance = case_when(
        p_val_adj < p_val_threshold & avg_log2FC > fc_threshold  ~ "Up in Cold",
        p_val_adj < p_val_threshold & avg_log2FC < -fc_threshold ~ "Down in Cold",
        TRUE ~ "Not Sig"
      )
    )
  
  # 处理 p_val_adj = 0，方便后续画图
  nonzero_pvals <- deg_results$p_val_adj[
    !is.na(deg_results$p_val_adj) & deg_results$p_val_adj > 0
  ]
  
  min_nonzero_pval <- ifelse(
    length(nonzero_pvals) > 0,
    min(nonzero_pvals),
    .Machine$double.xmin
  )
  
  deg_results$p_val_adj_plot <- deg_results$p_val_adj
  deg_results$p_val_adj_plot[deg_results$p_val_adj_plot == 0] <- min_nonzero_pval
  
  # 输出 CSV
  output_file <- file.path(
    "files",
    paste0("DEG_Aorta_", comparison_name, "_Cold_vs_RT.csv")
  )
  
  write.csv(
    deg_results,
    file = output_file,
    row.names = FALSE
  )
  
  print(paste0("✅ 已输出: ", output_file))
  
  return(deg_results)
}

# ------------------------------------------------------------------------------
# 3. 分别运行 4 组比较
# ------------------------------------------------------------------------------

deg_cluster7 <- run_cluster_deg(
  obj = aorta_obj,
  clusters_use = c("7"),
  comparison_name = "new_cluster7"
)

deg_cluster13 <- run_cluster_deg(
  obj = aorta_obj,
  clusters_use = c("13"),
  comparison_name = "new_cluster13"
)

deg_cluster16 <- run_cluster_deg(
  obj = aorta_obj,
  clusters_use = c("16"),
  comparison_name = "new_cluster16"
)

deg_cluster7_13 <- run_cluster_deg(
  obj = aorta_obj,
  clusters_use = c("7", "13"),
  comparison_name = "new_cluster7_13"
)
deg_cluster7_13_16 <- run_cluster_deg(
  obj = aorta_obj,
  clusters_use = c("7", "13", "16"),
  comparison_name = "new_cluster7_13_16"
)
print("🎉 全部 DEG 分析完成！")
# ==============================================================================
# 4. 对 cluster 7+13+16 的 DEG 结果进行 MSigDB Mouse GSEA
# ==============================================================================

library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)

if (!dir.exists("plots")) dir.create("plots")

print("🚀 开始 cluster 7+13+16 的 GSEA 富集分析...")

# ------------------------------------------------------------------------------
# 4.1 读取 MSigDB Mouse GO:BP 基因集
# ------------------------------------------------------------------------------

print("🚀 读取 MSigDB Mouse GO:BP gene sets...")

msigdb_gobp <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "GO:BP"
)

term2gene_gobp <- msigdb_gobp %>%
  dplyr::select(gs_name, gene_symbol)

# ------------------------------------------------------------------------------
# 4.2 构建 GSEA 排序基因列表
# ------------------------------------------------------------------------------

gsea_gene_list <- deg_cluster7_13_16 %>%
  dplyr::filter(!is.na(avg_log2FC)) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::select(gene, avg_log2FC)

gene_list <- gsea_gene_list$avg_log2FC
names(gene_list) <- gsea_gene_list$gene

gene_list <- sort(gene_list, decreasing = TRUE)

# ------------------------------------------------------------------------------
# 4.3 运行 GSEA
# ------------------------------------------------------------------------------

print("🚀 运行 GSEA...")

gsea_gobp <- GSEA(
  geneList = gene_list,
  TERM2GENE = term2gene_gobp,
  pvalueCutoff = 1,
  minGSSize = 10,
  maxGSSize = 500,
  eps = 0,
  verbose = FALSE
)

gsea_gobp_df <- as.data.frame(gsea_gobp)

# ------------------------------------------------------------------------------
# 4.4 保存 GSEA 结果表
# ------------------------------------------------------------------------------

write.csv(
  gsea_gobp_df,
  file = file.path("files", "GSEA_GO_BP_Aorta_new_cluster7_13_16_Cold_vs_RT.csv"),
  row.names = FALSE
)

print("✅ GSEA 结果表已保存。")

# ------------------------------------------------------------------------------
# 4.5 提取显著通路
# ------------------------------------------------------------------------------

gsea_sig <- gsea_gobp_df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(abs(NES)))

write.csv(
  gsea_sig,
  file = file.path("files", "GSEA_GO_BP_Aorta_new_cluster7_13_16_Cold_vs_RT_significant.csv"),
  row.names = FALSE
)

print(paste0("✅ 显著 GSEA 通路数量: ", nrow(gsea_sig)))

# ------------------------------------------------------------------------------
# 4.6 GSEA dotplot
# ------------------------------------------------------------------------------

if (nrow(gsea_gobp_df) > 0) {
  
  p_dot <- dotplot(
    gsea_gobp,
    showCategory = 20,
    split = ".sign"
  ) +
    facet_grid(. ~ .sign) +
    ggtitle("GSEA GO:BP - Aorta cluster 7+13+16 Cold vs RT")
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_dotplot_Aorta_new_cluster7_13_16_Cold_vs_RT.pdf"),
    plot = p_dot,
    width = 12,
    height = 8
  )
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_dotplot_Aorta_new_cluster7_13_16_Cold_vs_RT.png"),
    plot = p_dot,
    width = 12,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 4.7 NES barplot
# ------------------------------------------------------------------------------

top_gsea_plot <- gsea_gobp_df %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::arrange(desc(abs(NES))) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    Direction = ifelse(NES > 0, "Enriched in Cold", "Enriched in RT"),
    Description_short = stringr::str_replace_all(Description, "GOBP_", ""),
    Description_short = stringr::str_replace_all(Description_short, "_", " "),
    Description_short = stringr::str_to_title(Description_short),
    Description_short = stringr::str_wrap(Description_short, width = 45),
    Description_short = forcats::fct_reorder(Description_short, NES)
  )

if (nrow(top_gsea_plot) > 0) {
  
  p_bar <- ggplot(
    top_gsea_plot,
    aes(x = NES, y = Description_short, fill = Direction)
  ) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    theme_bw(base_size = 12) +
    labs(
      title = "Top GSEA GO:BP pathways",
      subtitle = "Aorta new cluster 7+13+16: Cold vs RT",
      x = "Normalized Enrichment Score, NES",
      y = NULL
    ) +
    theme(
      legend.title = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      axis.text.y = element_text(size = 10)
    )
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_NES_barplot_Aorta_new_cluster7_13_16_Cold_vs_RT.pdf"),
    plot = p_bar,
    width = 10,
    height = 8
  )
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_NES_barplot_Aorta_new_cluster7_13_16_Cold_vs_RT.png"),
    plot = p_bar,
    width = 10,
    height = 8,
    dpi = 300
  )
}

# ------------------------------------------------------------------------------
# 4.8 单个代表性通路 GSEA 曲线
# ------------------------------------------------------------------------------

if (nrow(gsea_sig) > 0) {
  
  top_pathway_id <- gsea_sig$ID[1]
  
  p_gseaplot <- gseaplot2(
    gsea_gobp,
    geneSetID = top_pathway_id,
    title = gsea_sig$Description[1],
    pvalue_table = TRUE
  )
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_top_pathway_Aorta_new_cluster7_13_16_Cold_vs_RT.pdf"),
    plot = p_gseaplot,
    width = 10,
    height = 7
  )
  
  ggsave(
    filename = file.path("plots", "GSEA_GO_BP_top_pathway_Aorta_new_cluster7_13_16_Cold_vs_RT.png"),
    plot = p_gseaplot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

print("🎉 cluster 7+13+16 GSEA 分析完成！")
# ==============================================================================
# Volcano Plot
# ==============================================================================

library(ggrepel)

print("🚀 绘制 Volcano Plot...")

# ------------------------------------------------------------------------------
# 选择需要标注的 top genes
# ------------------------------------------------------------------------------

top_up_genes <- deg_cluster7_13_16 %>%
  dplyr::filter(
    Significance == "Up in Cold"
  ) %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::slice_head(n = 10)

top_down_genes <- deg_cluster7_13_16 %>%
  dplyr::filter(
    Significance == "Down in Cold"
  ) %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::slice_head(n = 10)

top_genes <- bind_rows(
  top_up_genes,
  top_down_genes
)

# ------------------------------------------------------------------------------
# 火山图配色
# ------------------------------------------------------------------------------

volcano_colors <- c(
  "Up in Cold" = "#D51F26",
  "Down in Cold" = "#3C5488FF",
  "Not Sig" = "grey75"
)
deg_cluster7_13_16$Significance <- factor(
  deg_cluster7_13_16$Significance,
  levels = c(
    "Up in Cold",
    "Down in Cold",
    "Not Sig"
  )
)
# ------------------------------------------------------------------------------
# Volcano Plot
# ------------------------------------------------------------------------------

p_volcano <- ggplot(
  deg_cluster7_13_16,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj_plot),
    color = Significance
  )
) +
  
  geom_point(
    size = 1.2,
    alpha = 0.8
  ) +
  
  scale_color_manual(values = volcano_colors) +
  
  geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_hline(
    yintercept = -log10(p_val_threshold),
    linetype = "dashed",
    color = "grey40"
  ) +
  
  geom_text_repel(
    data = top_genes,
    aes(label = gene),
    size = 3.5,
    max.overlaps = Inf,
    box.padding = 0.5,
    point.padding = 0.3,
    show.legend = FALSE
  ) +
  
  labs(
    title = "Aorta new cluster 7+13+16",
    subtitle = "Cold vs RT",
    x = "avg_log2FC",
    y = "-log10(adjusted p value)"
  ) +
  
  theme_bw(base_size = 13) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    
    legend.title = element_blank()
  ) +
  
  coord_cartesian(
    xlim = c(-4, 4)
  )

# ------------------------------------------------------------------------------
# 保存图片
# ------------------------------------------------------------------------------

ggsave(
  filename = file.path(
    "plots",
    "Volcano_Aorta_new_cluster7_13_16_Cold_vs_RT.pdf"
  ),
  plot = p_volcano,
  width = 8,
  height = 7
)

ggsave(
  filename = file.path(
    "plots",
    "Volcano_Aorta_new_cluster7_13_16_Cold_vs_RT.png"
  ),
  plot = p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

print("✅ Volcano Plot 已输出。")