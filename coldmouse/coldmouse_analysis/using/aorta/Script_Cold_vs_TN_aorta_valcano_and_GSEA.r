# ==============================================================================
# Script: Whole Aorta Cold vs TN 差异分析 + MSigDB Mouse GSEA
# 重点输出图片
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggrepel)

if (!dir.exists("files")) dir.create("files")
if (!dir.exists("plots")) dir.create("plots")

# ------------------------------------------------------------------------------
# 0. 参数设置
# ------------------------------------------------------------------------------

tissue_name <- "Aorta"
analysis_scope <- "WholeAorta"

group_cold <- "Cold_4C"
group_tn   <- "TN_30C"

output_prefix <- paste(tissue_name, analysis_scope, "Cold_vs_TN", sep = "_")

fc_threshold <- 0.25
p_val_threshold <- 0.05

save_tables <- TRUE

print(paste0("📌 当前分析对象: ", tissue_name, " - ", analysis_scope))
print(paste0("📌 当前比较: ", group_cold, " vs ", group_tn))
print(paste0("📌 输出前缀: ", output_prefix))

# ------------------------------------------------------------------------------
# 1. 读取完整 Aorta 对象
# ------------------------------------------------------------------------------

print("🚀 步骤1: 读取完整 Aorta 对象...")

aorta_obj <- qread("aorta_corrected.qs")
DefaultAssay(aorta_obj) <- "RNA"

if (!"Group" %in% colnames(aorta_obj@meta.data)) {
  stop("metadata 中没有 Group 列，请检查对象。")
}

# 过滤 TN_30C，只保留 Cold_4C 和 TN_30C
aorta_obj <- subset(aorta_obj, subset = Group %in% c(group_cold, group_tn))

Idents(aorta_obj) <- "Group"

missing_groups <- setdiff(c(group_cold, group_tn), levels(Idents(aorta_obj)))

if (length(missing_groups) > 0) {
  stop(paste0("缺少以下分组: ", paste(missing_groups, collapse = ", ")))
}

print("  ✅ Aorta 对象准备完成。")
print(table(aorta_obj$Group))

# ------------------------------------------------------------------------------
# 2. 差异表达分析
# ------------------------------------------------------------------------------

print("🚀 步骤2: 计算 Whole Aorta Cold vs TN 差异基因...")

deg_results <- FindMarkers(
  aorta_obj,
  ident.1 = group_cold,
  ident.2 = group_tn,
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

# 先生成 p_val_adj_plot，再生成 top_genes，避免 geom_text_repel 报错
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

if (save_tables) {
  write.csv(
    deg_results,
    file = file.path("files", paste0("DEG_", output_prefix, ".csv")),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# 3. 火山图
# ------------------------------------------------------------------------------

print("🚀 步骤3: 绘制火山图...")

top_genes <- deg_results %>%
  filter(Significance != "Not Sig") %>%
  group_by(Significance) %>%
  slice_max(n = 10, order_by = abs(avg_log2FC)) %>%
  ungroup()

p_volcano <- ggplot(
  deg_results,
  aes(x = avg_log2FC, y = -log10(p_val_adj_plot), color = Significance)
) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(
    values = c(
      "Up in Cold" = "#D53E4F",
      "Down in Cold" = "#3288BD",
      "Not Sig" = "grey80"
    )
  ) +
  geom_vline(
    xintercept = c(-fc_threshold, fc_threshold),
    linetype = "dashed",
    color = "black"
  ) +
  geom_hline(
    yintercept = -log10(p_val_threshold),
    linetype = "dashed",
    color = "black"
  ) +
  geom_text_repel(
    data = top_genes,
    aes(
      x = avg_log2FC,
      y = -log10(p_val_adj_plot),
      label = gene
    ),
    inherit.aes = FALSE,
    size = 4,
    box.padding = 0.5,
    point.padding = 0.3,
    max.overlaps = 20,
    color = "black"
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = paste0("Volcano: ", tissue_name, " ", analysis_scope, " ", group_cold, " vs ", group_tn),
    x = expression("Log"[2] * " Fold Change"),
    y = expression("-Log"[10] * " Adjusted P-value")
  ) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

ggsave(
  filename = file.path("plots", paste0("Volcano_", output_prefix, ".png")),
  plot = p_volcano,
  width = 8,
  height = 7,
  dpi = 300
)

ggsave(
  filename = file.path("plots", paste0("Volcano_", output_prefix, ".pdf")),
  plot = p_volcano,
  width = 8,
  height = 7
)

# ------------------------------------------------------------------------------
# 4. 读取 MSigDB Mouse GO:BP 通路集
# ------------------------------------------------------------------------------

print("🚀 步骤4: 读取 MSigDB Mouse GO:BP 通路集...")

get_mouse_msigdb_gobp <- function() {
  
  msig_obj <- tryCatch(
    {
      msigdbr(
        db_species = "MM",
        species = "Mus musculus",
        collection = "M5",
        subcollection = "GO:BP"
      )
    },
    error = function(e1) {
      message("新版 msigdbr 参数失败，尝试旧版 M5 参数...")
      
      tryCatch(
        {
          msigdbr(
            species = "Mus musculus",
            category = "M5",
            subcategory = "GO:BP"
          )
        },
        error = function(e2) {
          message("旧版 M5 参数失败，尝试 C5 GO:BP 参数...")
          
          msigdbr(
            species = "Mus musculus",
            category = "C5",
            subcategory = "GO:BP"
          )
        }
      )
    }
  )
  
  return(as.data.frame(msig_obj))
}

msig_mouse_gobp <- get_mouse_msigdb_gobp()

if (nrow(msig_mouse_gobp) == 0) {
  stop("没有读取到 MSigDB mouse GO:BP 通路，请检查 msigdbr 包或网络环境。")
}

print("  当前 msigdbr 返回的列名如下：")
print(colnames(msig_mouse_gobp))

# 自动识别通路名列和基因名列
term_candidates <- c("gs_name", "gs_id", "gs_exact_source", "gs_standard_name", "gs_name_standard")
gene_candidates <- c("gene_symbol", "db_gene_symbol", "symbol", "mgi_symbol", "mouse_gene_symbol")

term_col <- intersect(term_candidates, colnames(msig_mouse_gobp))[1]
gene_col <- intersect(gene_candidates, colnames(msig_mouse_gobp))[1]

if (is.na(term_col) || is.na(gene_col)) {
  stop(
    paste0(
      "无法自动识别 MSigDB 的通路名列或基因名列。\n",
      "当前列名为：\n",
      paste(colnames(msig_mouse_gobp), collapse = ", ")
    )
  )
}

print(paste0("  ✅ 使用通路名列: ", term_col))
print(paste0("  ✅ 使用基因名列: ", gene_col))

term2gene <- msig_mouse_gobp %>%
  transmute(
    term = .data[[term_col]],
    gene = .data[[gene_col]]
  ) %>%
  distinct() %>%
  filter(!is.na(term), !is.na(gene))

# TERM2NAME 可有可无，但有的话图里显示更完整

term2name <- term2gene %>%
  distinct(term) %>%
  mutate(name = term)


print(paste0("  ✅ 读取到通路数量: ", length(unique(term2gene$term))))
print(paste0("  ✅ 基因-通路对应关系数量: ", nrow(term2gene)))

# ------------------------------------------------------------------------------
# 5. 构建 GSEA ranked gene list
# ------------------------------------------------------------------------------

print("🚀 步骤5: 构建 GSEA ranked gene list...")

gene_list_df <- deg_results %>%
  dplyr::select(gene, avg_log2FC) %>%
  filter(!is.na(gene), !is.na(avg_log2FC)) %>%
  group_by(gene) %>%
  summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop") %>%
  arrange(desc(avg_log2FC))

gene_list_vector <- gene_list_df$avg_log2FC
names(gene_list_vector) <- gene_list_df$gene
gene_list_vector <- sort(gene_list_vector, decreasing = TRUE)

# 只保留 MSigDB mouse 中存在的基因
gene_list_vector <- gene_list_vector[
  names(gene_list_vector) %in% unique(term2gene$gene)
]

print(paste0("  ✅ 进入 GSEA 的基因数: ", length(gene_list_vector)))

if (length(gene_list_vector) < 100) {
  stop("进入 GSEA 的基因数少于 100，请检查基因名是否为小鼠 SYMBOL。")
}

# ------------------------------------------------------------------------------
# 6. 运行 MSigDB Mouse GSEA
# ------------------------------------------------------------------------------

print("🚀 步骤6: 运行 MSigDB Mouse GSEA...")

set.seed(123)

gsea_res <- GSEA(
  geneList      = gene_list_vector,
  TERM2GENE     = term2gene,
  TERM2NAME     = term2name,
  minGSSize     = 10,
  maxGSSize     = 500,
  pvalueCutoff  = 1,
  pAdjustMethod = "BH",
  eps           = 0,
  verbose       = FALSE,
  seed          = TRUE,
  by            = "fgsea"
)

gsea_df <- as.data.frame(gsea_res)

if (save_tables) {
  write.csv(
    gsea_df,
    file = file.path("files", paste0("GSEA_MSigDB_Mouse_GOBP_", output_prefix, ".csv")),
    row.names = FALSE
  )
}

print(paste0("  ✅ GSEA 通路总数: ", nrow(gsea_df)))
print(paste0("  ✅ 显著通路数量 p.adjust < 0.05: ", sum(gsea_df$p.adjust < 0.05, na.rm = TRUE)))

# ------------------------------------------------------------------------------
# 7. 美化版 GSEA dotplot
# ------------------------------------------------------------------------------

print("🚀 步骤7: 绘制美化版 GSEA dotplot...")

# 美化通路名称函数
beautify_pathway <- function(x) {
  x %>%
    stringr::str_remove("^GOBP_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title() %>%
    stringr::str_wrap(width = 32)
}

# 从 GSEA 结果中提取 top pathways
# NES > 0: Cold enriched
# NES < 0: TN enriched
top_gsea_dot <- bind_rows(
  gsea_df %>%
    filter(NES > 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10),
  
  gsea_df %>%
    filter(NES < 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
) %>%
  mutate(
    Direction = ifelse(
      NES > 0,
      "Enriched in Cold_4C",
      "Enriched in TN_30C"
    ),
    
    # 用 ID 而不是 Description，避免超长 GO 描述
    Pathway = beautify_pathway(ID),
    
    # 如果 Count 不存在，则用 core_enrichment 计算 leading-edge gene 数
    Count = ifelse(
      "Count" %in% colnames(.),
      Count,
      lengths(strsplit(core_enrichment, "/"))
    ),
    
    GeneRatio = Count / setSize
  )

# 按 NES 排序，让 Cold 和 TN 两侧更有方向感
top_gsea_dot <- top_gsea_dot %>%
  arrange(NES) %>%
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)))

p_dot <- ggplot(
  top_gsea_dot,
  aes(
    x = GeneRatio,
    y = Pathway,
    size = Count,
    color = p.adjust
  )
) +
  geom_point(alpha = 0.9) +
  facet_grid(. ~ Direction) +
  scale_color_gradient(
    low = "#D53E4F",
    high = "#3288BD",
    name = "Adjusted P"
  ) +
  scale_size_continuous(name = "Gene Count") +
  theme_bw(base_size = 13) +
  labs(
    title = paste0("MSigDB Mouse GO:BP GSEA: ", output_prefix),
    x = "Gene Ratio",
    y = NULL
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path("plots", paste0("GSEA_Dotplot_Beautified_", output_prefix, ".png")),
  plot = p_dot,
  width = 12,
  height = 10,
  dpi = 300
)

ggsave(
  filename = file.path("plots", paste0("GSEA_Dotplot_Beautified_", output_prefix, ".pdf")),
  plot = p_dot,
  width = 12,
  height = 10
)

print("  ✅ 美化版 GSEA dotplot 已保存。")

# ------------------------------------------------------------------------------
# 8. GSEA NES barplot
# ------------------------------------------------------------------------------

print("🚀 步骤8: 绘制 GSEA NES barplot...")

top_gsea_bar <- bind_rows(
  gsea_df %>%
    filter(NES > 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10),
  gsea_df %>%
    filter(NES < 0) %>%
    arrange(p.adjust) %>%
    slice_head(n = 10)
) %>%
  mutate(
    Direction = ifelse(NES > 0, "Enriched in Cold_4C", "Enriched in TN_30C"),
    
    # 使用 ID，而不是 Description，避免超长 GO 描述
    Pathway = ID,
    
    # 美化通路名
    Pathway = str_remove(Pathway, "^GOBP_"),
    Pathway = str_replace_all(Pathway, "_", " "),
    Pathway = str_to_title(Pathway),
    
    # 自动换行，防止标签过长
    Pathway = str_wrap(Pathway, width = 35),
    
    # 按 NES 排序
    Pathway = factor(Pathway, levels = unique(Pathway[order(NES)]))
  )

if (nrow(top_gsea_bar) > 0) {
  
  p_bar <- ggplot(
    top_gsea_bar,
    aes(x = NES, y = Pathway, fill = Direction)
  ) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    scale_fill_manual(
      values = c(
        "Enriched in Cold_4C" = "#D53E4F",
        "Enriched in TN_30C" = "#3288BD"
      )
    ) +
    theme_classic(base_size = 13) +
    labs(
      title = paste0("Top MSigDB Mouse GO:BP GSEA Terms: ", output_prefix),
      x = "Normalized Enrichment Score NES",
      y = NULL
    ) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.y = element_text(size = 9)
    )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_NES_Barplot_", output_prefix, ".png")),
    plot = p_bar,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_NES_Barplot_", output_prefix, ".pdf")),
    plot = p_bar,
    width = 12,
    height = 8
  )
}
# ------------------------------------------------------------------------------
# 9. Top Cold / TN 富集通路曲线图
# ------------------------------------------------------------------------------

print("🚀 步骤9: 绘制 top GSEA enrichment curve...")

top_cold_term <- gsea_df %>%
  filter(NES > 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 1)

top_rt_term <- gsea_df %>%
  filter(NES < 0) %>%
  arrange(p.adjust) %>%
  slice_head(n = 1)

if (nrow(top_cold_term) == 1) {
  
  p_top_cold <- gseaplot2(
    gsea_res,
    geneSetID = top_cold_term$ID[1],
    title = paste0("Top Cold-enriched: ", top_cold_term$Description[1])
  )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_Top_Cold_", output_prefix, ".png")),
    plot = p_top_cold,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_Top_Cold_", output_prefix, ".pdf")),
    plot = p_top_cold,
    width = 8,
    height = 6
  )
}

if (nrow(top_rt_term) == 1) {
  
  p_top_rt <- gseaplot2(
    gsea_res,
    geneSetID = top_rt_term$ID[1],
    title = paste0("Top TN-enriched: ", top_rt_term$Description[1])
  )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_Top_TN_", output_prefix, ".png")),
    plot = p_top_rt,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path("plots", paste0("GSEA_Top_TN_", output_prefix, ".pdf")),
    plot = p_top_rt,
    width = 8,
    height = 6
  )
}

print("🎉 Whole Aorta GSEA 分析完成！")
print("重点图片输出:")
print(paste0("  plots/Volcano_", output_prefix, ".png"))
print(paste0("  plots/GSEA_Dotplot_", output_prefix, ".png"))
print(paste0("  plots/GSEA_NES_Barplot_", output_prefix, ".png"))
print(paste0("  plots/GSEA_Top_Cold_", output_prefix, ".png"))
print(paste0("  plots/GSEA_Top_TN_", output_prefix, ".png"))