# ==============================================================================
# 05_marker_and_enrichment_introns.R
# coldmouse_introns：DEG + Volcano + GO:BP GSEA 分析
# ==============================================================================
# 输入：coldmouse_introns_sc_by_tissue_no_harmony_louvain.qs
# 输出：DEG csv、Volcano png/pdf、GSEA csv/png/pdf
# 参考：Script_Cold_vs_RT_aorta_valcano_and_GSEA.r / Script_Cold_vs_TN_aorta_valcano_and_GSEA.r
# ==============================================================================

PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))

config <- load_config()
setup_performance(config)

library(Seurat)
library(tidyverse)
library(qs)
library(ggplot2)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(patchwork)

# 输出目录
output_dir <- config$paths$output_dir
files_dir <- file.path(output_dir, "files")
plots_dir <- file.path(output_dir, "plots")
dir.create(files_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

# 参数
fc_threshold <- 0.25
p_val_threshold <- 0.05

# 读取按组织对象
sc_by_tissue_file <- file.path(output_dir,
                               paste0("coldmouse_introns_sc_by_tissue_no_harmony_",
                                      config$clustering$method, ".qs"))
sc_by_tissue <- qread(sc_by_tissue_file)

# ------------------------------------------------------------------------------
# 辅助函数
# ------------------------------------------------------------------------------

# 获取 MSigDB GO:BP 基因集 (三级回退)
get_gobp_genesets <- function() {
  cat("获取 MSigDB GO:BP 基因集...\n")

  term2gene <- NULL
  term2name <- NULL

  # 尝试1: collection + subcollection (新 API)
  tryCatch({
    gs <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
    if (nrow(gs) > 0) {
      term_col <- intersect(c("gs_name", "gene_symbol"), colnames(gs))
      term2gene <- gs %>% select(all_of(term_col))
      colnames(term2gene) <- c("term", "gene")
      cat("  使用 collection='C5', subcollection='GO:BP'\n")
    }
  }, error = function(e) cat("  尝试1失败:", e$message, "\n"))

  # 尝试2: category + subcategory (旧 API)
  if (is.null(term2gene)) {
    tryCatch({
      gs <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP")
      if (nrow(gs) > 0) {
        term_col <- intersect(c("gs_name", "gene_symbol"), colnames(gs))
        term2gene <- gs %>% select(all_of(term_col))
        colnames(term2gene) <- c("term", "gene")
        cat("  使用 category='C5', subcategory='GO:BP'\n")
      }
    }, error = function(e) cat("  尝试2失败:", e$message, "\n"))
  }

  if (is.null(term2gene)) {
    stop("无法获取 MSigDB GO:BP 基因集")
  }

  # term2name
  term2name <- term2gene %>%
    distinct(term) %>%
    mutate(name = gsub("^GOBP_", "", term) %>%
             gsub("_", " ", .) %>%
             tools::toTitleCase())

  list(term2gene = term2gene, term2name = term2name)
}

# 火山图绘制
plot_volcano <- function(deg_df, title, fc_thresh = 0.25, p_thresh = 0.05,
                         top_n = 10, width = 10, height = 8, dpi = 300) {
  # 处理零 p 值
  min_nonzero_pval <- min(deg_df$p_val[deg_df$p_val > 0], na.rm = TRUE)
  deg_df$p_val_adj_safe <- pmax(deg_df$p_val_adj, 1e-300)
  deg_df$p_val_safe <- pmax(deg_df$p_val, min_nonzero_pval)

  deg_df$direction <- case_when(
    deg_df$avg_log2FC > fc_thresh & deg_df$p_val < p_thresh ~ "Up in Cold",
    deg_df$avg_log2FC < -fc_thresh & deg_df$p_val < p_thresh ~ "Down in Cold",
    TRUE ~ "Not Sig"
  )

  # 标记 top 基因
  up_genes <- deg_df %>% filter(direction == "Up in Cold") %>%
    arrange(desc(avg_log2FC)) %>% head(top_n)
  down_genes <- deg_df %>% filter(direction == "Down in Cold") %>%
    arrange(avg_log2FC) %>% head(top_n)
  label_genes <- bind_rows(up_genes, down_genes) %>% distinct(gene)

  deg_df$label <- ifelse(deg_df$gene %in% label_genes$gene, deg_df$gene, "")

  p <- ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_safe), color = direction)) +
    geom_point(size = 0.8, alpha = 0.6) +
    scale_color_manual(values = c("Up in Cold" = "#D53E4F", "Down in Cold" = "#3288BD", "Not Sig" = "grey80")) +
    geom_vline(xintercept = c(-fc_thresh, fc_thresh), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(p_thresh), linetype = "dashed", color = "grey40") +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 20, show.legend = FALSE) +
    ggtitle(title) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")

  list(plot = p, data = deg_df)
}

# GSEA 点图
plot_gsea_dotplot <- function(gsea_res, title, top_n = 10) {
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    cat("  GSEA 结果为空，跳过绘图\n")
    return(NULL)
  }

  # 美化 pathway 名称
  clean_name <- function(x) {
    x <- gsub("^GOBP_", "", x)
    x <- gsub("_", " ", x)
    x <- tools::toTitleCase(x)
    x
  }

  gsea_res$Description_clean <- clean_name(gsea_res$Description)

  # 分方向取 top
  up_pathways <- gsea_res %>% filter(NES > 0) %>%
    arrange(p.adjust) %>% head(top_n)
  down_pathways <- gsea_res %>% filter(NES < 0) %>%
    arrange(p.adjust) %>% head(top_n)
  plot_data <- bind_rows(up_pathways, down_pathways) %>%
    mutate(Direction = ifelse(NES > 0, "Enriched in Cold", "Enriched in Control"))

  plot_data$GeneRatio <- plot_data$setSize / max(plot_data$setSize)

  p <- ggplot(plot_data, aes(x = GeneRatio, y = reorder(Description_clean, GeneRatio))) +
    geom_point(aes(color = p.adjust, size = setSize)) +
    scale_color_gradient(low = "#D53E4F", high = "#3288BD") +
    facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5, face = "bold"))

  p
}

# NES 柱状图
plot_gsea_barplot <- function(gsea_res, title, top_n = 10) {
  if (is.null(gsea_res) || nrow(gsea_res) == 0) return(NULL)

  clean_name <- function(x) {
    x <- gsub("^GOBP_", "", x)
    x <- gsub("_", " ", x)
    tools::toTitleCase(x)
  }

  up_pathways <- gsea_res %>% filter(NES > 0) %>%
    arrange(p.adjust) %>% head(top_n)
  down_pathways <- gsea_res %>% filter(NES < 0) %>%
    arrange(p.adjust) %>% head(top_n)
  plot_data <- bind_rows(up_pathways, down_pathways) %>%
    mutate(Direction = ifelse(NES > 0, "Enriched in Cold", "Enriched in Control"),
           Description_clean = clean_name(Description))

  ggplot(plot_data, aes(x = NES, y = reorder(Description_clean, NES), fill = Direction)) +
    geom_col() +
    scale_fill_manual(values = c("Enriched in Cold" = "#D53E4F", "Enriched in Control" = "#3288BD")) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 9),
          plot.title = element_text(hjust = 0.5, face = "bold"))
}

# 单个比较分析
run_deg_gsea <- function(obj, tissue_name, group1, group2, output_prefix,
                         gobp_genesets, files_dir, plots_dir) {
  cat("\n===================================================================\n")
  cat("分析:", output_prefix, "\n")
  cat("===================================================================\n")

  # 子集
  obj_sub <- subset(obj, subset = Group %in% c(group1, group2))
  Idents(obj_sub) <- "Group"

  cat("  细胞数:", ncol(obj_sub), "(", group1, ":",
      sum(obj_sub$Group == group1), ",", group2, ":",
      sum(obj_sub$Group == group2), ")\n")

  # DEG
  cat("  计算 DEG...\n")
  deg <- FindMarkers(obj_sub, ident.1 = group1, ident.2 = group2,
                     logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  deg$gene <- rownames(deg)

  # 保存 DEG
  write.csv(deg, file.path(files_dir, paste0("introns_DEG_", output_prefix, ".csv")),
            row.names = FALSE)
  cat("  DEG 表已保存 (", nrow(deg), "个基因 )\n")

  # 火山图
  volcano_result <- plot_volcano(deg, paste0("Introns ", output_prefix))
  ggsave(file.path(plots_dir, paste0("introns_Volcano_", output_prefix, ".png")),
         plot = volcano_result$plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(plots_dir, paste0("introns_Volcano_", output_prefix, ".pdf")),
         plot = volcano_result$plot, width = 10, height = 8)

  # 统计显著基因
  sig_up <- sum(deg$avg_log2FC > fc_threshold & deg$p_val < p_val_threshold)
  sig_down <- sum(deg$avg_log2FC < -fc_threshold & deg$p_val < p_val_threshold)
  cat("  显著上调基因:", sig_up, "显著下调基因:", sig_down, "\n")

  # GSEA
  cat("  运行 GO:BP GSEA...\n")
  ranked_list <- deg %>%
    select(gene, avg_log2FC) %>%
    group_by(gene) %>%
    summarise(avg_log2FC = mean(avg_log2FC), .groups = "drop") %>%
    arrange(desc(avg_log2FC))

  # 过滤到基因集中的基因
  genes_in_sets <- unique(gobp_genesets$term2gene$gene)
  ranked_list <- ranked_list %>% filter(gene %in% genes_in_sets)

  gene_vec <- setNames(ranked_list$avg_log2FC, ranked_list$gene)

  gsea_res <- NULL
  tryCatch({
    gsea_res <- GSEA(geneList = gene_vec,
                     TERM2GENE = gobp_genesets$term2gene,
                     TERM2NAME = gobp_genesets$term2name,
                     by = "fgsea",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 1,
                     pAdjustMethod = "BH",
                     eps = 0,
                     verbose = FALSE)

    if (!is.null(gsea_res) && nrow(as.data.frame(gsea_res)) > 0) {
      gsea_df <- as.data.frame(gsea_res)
      write.csv(gsea_df, file.path(files_dir, paste0("introns_GSEA_GOBP_", output_prefix, ".csv")),
                row.names = FALSE)
      cat("  GSEA 完成 (", nrow(gsea_df), "条通路 )\n")

      # GSEA 点图
      p_dot <- plot_gsea_dotplot(gsea_df, paste0("Introns GO:BP GSEA ", output_prefix))
      if (!is.null(p_dot)) {
        ggsave(file.path(plots_dir, paste0("introns_GSEA_Dotplot_", output_prefix, ".png")),
               plot = p_dot, width = 12, height = 10, dpi = 300)
        ggsave(file.path(plots_dir, paste0("introns_GSEA_Dotplot_", output_prefix, ".pdf")),
               plot = p_dot, width = 12, height = 10)
      }

      # NES 柱状图
      p_bar <- plot_gsea_barplot(gsea_df, paste0("Introns NES ", output_prefix))
      if (!is.null(p_bar)) {
        ggsave(file.path(plots_dir, paste0("introns_GSEA_NES_Barplot_", output_prefix, ".png")),
               plot = p_bar, width = 10, height = 8, dpi = 300)
      }

      # 富集曲线 (top Cold 和 top Control)
      tryCatch({
        up_terms <- gsea_df %>% filter(NES > 0) %>% arrange(p.adjust)
        down_terms <- gsea_df %>% filter(NES < 0) %>% arrange(p.adjust)

        if (nrow(up_terms) > 0) {
          p_up <- enrichplot::gseaplot2(gsea_res, geneSetID = up_terms$ID[1:3],
                            title = "Top Cold-enriched pathways")
          ggsave(file.path(plots_dir, paste0("introns_GSEA_Curve_Cold_", output_prefix, ".png")),
                 plot = p_up, width = 8, height = 8, dpi = 300)
        }
        if (nrow(down_terms) > 0) {
          p_down <- enrichplot::gseaplot2(gsea_res, geneSetID = down_terms$ID[1:3],
                              title = "Top Control-enriched pathways")
          ggsave(file.path(plots_dir, paste0("introns_GSEA_Curve_Control_", output_prefix, ".png")),
                 plot = p_down, width = 8, height = 8, dpi = 300)
        }
      }, error = function(e) cat("  富集曲线绘制失败:", e$message, "\n"))
    } else {
      cat("  GSEA 无显著结果\n")
    }
  }, error = function(e) {
    cat("  GSEA 失败:", e$message, "\n")
  })

  cat("  完成:", output_prefix, "\n")
}

# ------------------------------------------------------------------------------
# 主流程
# ------------------------------------------------------------------------------
cat("=== coldmouse_introns DEG + GSEA 分析 ===\n")

gobp_genesets <- get_gobp_genesets()
cat("GO:BP 基因集数量:", nrow(gobp_genesets$term2name), "\n")

# 对每个组织运行 Cold vs RT 和 Cold vs TN
for (tissue in names(sc_by_tissue)) {
  obj <- sc_by_tissue[[tissue]]

  # Cold vs RT
  run_deg_gsea(obj, tissue, "Cold_4C", "RT_25C",
               paste0(tissue, "_Cold_vs_RT"), gobp_genesets, files_dir, plots_dir)

  # Cold vs TN
  run_deg_gsea(obj, tissue, "Cold_4C", "TN_30C",
               paste0(tissue, "_Cold_vs_TN"), gobp_genesets, files_dir, plots_dir)
}

cat("\n=== 所有分析完成 ===\n")
cat("结果文件位于:", files_dir, "\n")
cat("图片文件位于:", plots_dir, "\n")
