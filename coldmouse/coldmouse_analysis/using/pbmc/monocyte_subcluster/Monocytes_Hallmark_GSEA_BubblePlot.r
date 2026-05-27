# ==============================================================================
# Script 5: Monocytes 亚群 Hallmark GSEA 分析 + 气泡图绘制
# ==============================================================================

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(qs)
  library(msigdbr)
  library(clusterProfiler)
})

set.seed(42)

# ------------------------------------------------------------------------------
# 0. 基础设置
# ------------------------------------------------------------------------------

dir.create("files", showWarnings = FALSE, recursive = TRUE)
dir.create("pictures", showWarnings = FALSE, recursive = TRUE)

mono_qs_file <- "pbmc_monocytes_sub-clustered.qs"

# 物种：你的项目是 mouse PBMC，这里使用 Mus musculus
species_name <- "Mus musculus"

# GSEA 参数
min_gs_size <- 10
max_gs_size <- 500
fdr_cutoff <- 0.25

# 气泡大小使用 -log10(FDR)
# 颜色使用 NES：红色代表该 mono cluster 中富集，上调；蓝色代表下调
plot_size_by <- "neg_log10_padj"


# ------------------------------------------------------------------------------
# 1. 读取 Monocytes 子集对象
# ------------------------------------------------------------------------------

print("🚀 步骤1: 读取 pbmc_monocytes_sub-clustered.qs ...")

if (!file.exists(mono_qs_file)) {
  stop(paste0("未找到文件: ", mono_qs_file))
}

mono_cells <- qread(mono_qs_file)

required_meta <- c("mono_cluster_id", "mono_annotation")
missing_meta <- setdiff(required_meta, colnames(mono_cells@meta.data))

if (length(missing_meta) > 0) {
  stop(paste0(
    "mono_cells@meta.data 中缺少以下列: ",
    paste(missing_meta, collapse = ", ")
  ))
}

# 优先使用 RNA assay 做差异分析
assay_to_use <- if ("RNA" %in% Assays(mono_cells)) "RNA" else DefaultAssay(mono_cells)
DefaultAssay(mono_cells) <- assay_to_use
print(paste0("  ✅ 使用 assay: ", assay_to_use))

# Seurat v5 多 layer 对象可能需要 JoinLayers
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  mono_cells <- tryCatch(
    {
      JoinLayers(mono_cells)
    },
    error = function(e) {
      message("  ⚠️ JoinLayers 跳过: ", e$message)
      mono_cells
    }
  )
}


# ------------------------------------------------------------------------------
# 2. 整理 Mono_1 到 Mono_6 的顺序和注释
# ------------------------------------------------------------------------------

print("🚀 步骤2: 整理 Mono_1 到 Mono_6 的注释信息...")

sort_mono_id <- function(x) {
  num <- suppressWarnings(as.numeric(stringr::str_extract(x, "\\d+")))
  x[order(num, x, na.last = TRUE)]
}

available_clusters <- sort_mono_id(unique(as.character(mono_cells$mono_cluster_id)))
target_clusters <- paste0("Mono_", 1:6)

missing_clusters <- setdiff(target_clusters, available_clusters)
if (length(missing_clusters) > 0) {
  stop(paste0(
    "对象中缺少以下 mono_cluster_id: ",
    paste(missing_clusters, collapse = ", ")
  ))
}

mono_cells$mono_cluster_id <- factor(
  as.character(mono_cells$mono_cluster_id),
  levels = target_clusters
)

Idents(mono_cells) <- "mono_cluster_id"

annotation_meta <- mono_cells@meta.data %>%
  as_tibble(rownames = "cell") %>%
  select(mono_cluster_id, mono_annotation) %>%
  mutate(
    mono_cluster_id = as.character(mono_cluster_id),
    mono_annotation = as.character(mono_annotation)
  ) %>%
  filter(mono_cluster_id %in% target_clusters) %>%
  group_by(mono_cluster_id) %>%
  summarise(
    mono_annotation = {
      x <- unique(mono_annotation[!is.na(mono_annotation) & mono_annotation != ""])
      if (length(x) == 0) unique(mono_cluster_id)[1] else x[1]
    },
    .groups = "drop"
  )

annotation_df <- tibble(mono_cluster_id = target_clusters) %>%
  left_join(annotation_meta, by = "mono_cluster_id") %>%
  mutate(
    mono_annotation = if_else(is.na(mono_annotation), mono_cluster_id, mono_annotation),
    mono_y_label =  mono_annotation
  )

write.csv(
  annotation_df,
  file.path("files", "Monocytes_MonoCluster_Annotation_For_GSEA.csv"),
  row.names = FALSE
)

print("  ✅ Mono cluster 注释如下:")
print(annotation_df)


# ------------------------------------------------------------------------------
# 3. 获取 MSigDB Hallmark 基因集
# ------------------------------------------------------------------------------

print("🚀 步骤3: 获取 MSigDB Hallmark gene sets ...")

get_hallmark_msigdb <- function(species_name) {
  
  msig_obj <- tryCatch(
    {
      msigdbr(species = species_name, category = "H")
    },
    error = function(e1) {
      msigdbr(species = species_name, collection = "H")
    }
  )
  
  term_col <- intersect(c("gs_name", "gs_id"), colnames(msig_obj))[1]
  gene_col <- intersect(c("gene_symbol", "db_gene_symbol", "human_gene_symbol"), colnames(msig_obj))[1]
  
  if (is.na(term_col) || is.na(gene_col)) {
    stop("msigdbr 返回结果中未找到 gs_name / gene_symbol 相关列，请检查 msigdbr 版本。")
  }
  
  out <- msig_obj[, c(term_col, gene_col)]
  colnames(out) <- c("gs_name", "gene_symbol")
  
  out %>%
    as_tibble() %>%
    filter(!is.na(gs_name), !is.na(gene_symbol)) %>%
    distinct(gs_name, gene_symbol)
}

term2gene <- get_hallmark_msigdb(species_name)

# 检查基因名大小写是否匹配
obj_genes <- rownames(mono_cells)

overlap_raw <- length(intersect(obj_genes, term2gene$gene_symbol))
overlap_upper <- length(intersect(toupper(obj_genes), toupper(term2gene$gene_symbol)))

use_uppercase <- FALSE

if (overlap_upper > overlap_raw * 2 && overlap_raw < 500) {
  print("  ⚠️ 检测到基因名大小写可能不匹配，将使用大写 gene symbol 进行 GSEA。")
  use_uppercase <- TRUE
  term2gene <- term2gene %>%
    mutate(gene_symbol = toupper(gene_symbol)) %>%
    distinct(gs_name, gene_symbol)
}

final_overlap <- if (use_uppercase) {
  length(intersect(toupper(obj_genes), term2gene$gene_symbol))
} else {
  length(intersect(obj_genes, term2gene$gene_symbol))
}

print(paste0("  ✅ Hallmark gene sets 数量: ", length(unique(term2gene$gs_name))))
print(paste0("  ✅ 与对象表达矩阵重叠的基因数: ", final_overlap))

if (final_overlap < 100) {
  stop("Hallmark 基因集与表达矩阵基因名重叠过少，请检查物种或 gene symbol 格式。")
}

write.csv(
  term2gene,
  file.path("files", "MSigDB_Hallmark_TERM2GENE_Used_For_GSEA.csv"),
  row.names = FALSE
)


# ------------------------------------------------------------------------------
# 4. 工具函数
# ------------------------------------------------------------------------------

get_logfc_col <- function(marker_df) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "avg_log2fc")
  fc_col <- intersect(fc_candidates, colnames(marker_df))[1]
  
  if (is.na(fc_col)) {
    stop("FindMarkers 结果中未找到 avg_log2FC 或 avg_logFC 列。")
  }
  
  fc_col
}

clean_hallmark_name <- function(x) {
  
  y <- x %>%
    stringr::str_remove("^HALLMARK_") %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_title()
  
  y <- stringr::str_replace_all(
    y,
    c(
      "Mtorc1" = "MTORC1",
      "Mtor" = "MTOR",
      "Pi3k" = "PI3K",
      "Akt" = "AKT",
      "Myc" = "MYC",
      "Kras" = "KRAS",
      "Il6" = "IL6",
      "Jak" = "JAK",
      "Stat3" = "STAT3",
      "Tnfa" = "TNFA",
      "Nfkb" = "NFKB",
      "Tgf" = "TGF",
      "Uv" = "UV",
      "Dna" = "DNA",
      "G2m" = "G2M",
      "P53" = "P53",
      "Ros" = "ROS"
    )
  )
  
  y
}


# ------------------------------------------------------------------------------
# 5. 对每个 Mono cluster 做 cluster vs others 的 ranked GSEA
# ------------------------------------------------------------------------------

print("🚀 步骤4: 开始对 Mono_1 到 Mono_6 分别进行 Hallmark GSEA ...")

run_one_cluster_gsea <- function(cluster_id) {
  
  print(paste0("  🔍 正在处理 ", cluster_id, " vs others ..."))
  
  # GSEA 不建议先过滤基因，所以 logfc.threshold = 0, min.pct = 0
  marker_df <- FindMarkers(
    mono_cells,
    ident.1 = cluster_id,
    ident.2 = NULL,
    only.pos = FALSE,
    min.pct = 0,
    logfc.threshold = 0,
    test.use = "wilcox",
    verbose = FALSE
  ) %>%
    rownames_to_column(var = "gene")
  
  fc_col <- get_logfc_col(marker_df)
  
  # 这里用 avg_log2FC / avg_logFC 作为排序统计量
  # 正值：该 cluster 相对其他 Monocytes 上调
  # 负值：该 cluster 相对其他 Monocytes 下调
  ranked_df <- marker_df %>%
    transmute(
      gene = if (use_uppercase) toupper(gene) else gene,
      rank_score = .data[[fc_col]]
    ) %>%
    filter(!is.na(gene), !is.na(rank_score)) %>%
    group_by(gene) %>%
    summarise(rank_score = mean(rank_score), .groups = "drop") %>%
    arrange(desc(rank_score), gene)
  
  gene_list <- ranked_df$rank_score
  names(gene_list) <- ranked_df$gene
  gene_list <- sort(gene_list, decreasing = TRUE)
  
  ranked_file <- file.path(
    "files",
    paste0("RankedGenes_", cluster_id, "_vs_Others_For_Hallmark_GSEA.csv")
  )
  
  write.csv(ranked_df, ranked_file, row.names = FALSE)
  
  gsea_obj <- suppressWarnings(
    GSEA(
      geneList = gene_list,
      TERM2GENE = term2gene,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      pvalueCutoff = 1,
      pAdjustMethod = "BH",
      eps = 0,
      verbose = FALSE,
      seed = TRUE
    )
  )
  
  gsea_df <- as.data.frame(gsea_obj) %>%
    as_tibble()
  
  if (nrow(gsea_df) == 0) {
    warning(paste0(cluster_id, " 没有得到 GSEA 结果。"))
    return(tibble())
  }
  
  gsea_df <- gsea_df %>%
    mutate(
      mono_cluster_id = cluster_id,
      pathway = ID,
      pathway_clean = clean_hallmark_name(ID),
      core_gene_count = if_else(
        is.na(core_enrichment) | core_enrichment == "",
        0L,
        lengths(strsplit(core_enrichment, "/", fixed = TRUE))
      ),
      core_gene_fraction = core_gene_count / setSize,
      p.adjust = if_else(is.na(p.adjust), 1, p.adjust),
      neg_log10_padj = pmin(-log10(p.adjust + 1e-300), 20),
      significant_FDR_0.25 = p.adjust <= fdr_cutoff
    ) %>%
    left_join(annotation_df, by = "mono_cluster_id") %>%
    relocate(
      mono_cluster_id,
      mono_annotation,
      mono_y_label,
      pathway,
      pathway_clean,
      NES,
      pvalue,
      p.adjust,
      qvalue,
      setSize,
      core_gene_count,
      core_gene_fraction,
      neg_log10_padj,
      significant_FDR_0.25
    )
  
  out_file <- file.path(
    "files",
    paste0("GSEA_Hallmark_", cluster_id, "_vs_Others.csv")
  )
  
  write.csv(gsea_df, out_file, row.names = FALSE)
  
  print(paste0("    ✅ 完成: ", cluster_id))
  
  gsea_df
}

gsea_combined <- map_dfr(target_clusters, run_one_cluster_gsea)

if (nrow(gsea_combined) == 0) {
  stop("所有 cluster 均未得到 GSEA 结果。")
}

combined_file <- file.path(
  "files",
  "GSEA_Hallmark_Mono_1_to_Mono_6_vs_Others_Combined.csv"
)

write.csv(gsea_combined, combined_file, row.names = FALSE)

sig_file <- file.path(
  "files",
  "GSEA_Hallmark_Mono_1_to_Mono_6_vs_Others_FDR0.25.csv"
)

write.csv(
  gsea_combined %>% filter(significant_FDR_0.25),
  sig_file,
  row.names = FALSE
)

print(paste0("  ✅ 合并 GSEA 结果已导出至: ", combined_file))


# ------------------------------------------------------------------------------
# 6. 整理气泡图数据
# ------------------------------------------------------------------------------

print("🚀 步骤5: 整理气泡图数据...")

plot_df <- gsea_combined %>%
  mutate(
    plot_size = .data[[plot_size_by]]
  )

# 按不同 Mono cluster 中 NES 的模式对 pathway 排序，便于观察相似通路模块
pathway_nes_mat_df <- plot_df %>%
  select(pathway_clean, mono_cluster_id, NES) %>%
  distinct() %>%
  tidyr::pivot_wider(
    names_from = mono_cluster_id,
    values_from = NES,
    values_fill = 0
  )

if (nrow(pathway_nes_mat_df) > 1) {
  pathway_nes_mat <- pathway_nes_mat_df %>%
    column_to_rownames("pathway_clean") %>%
    as.matrix()
  
  pathway_order <- rownames(pathway_nes_mat)[
    hclust(dist(pathway_nes_mat))$order
  ]
} else {
  pathway_order <- pathway_nes_mat_df$pathway_clean
}

plot_df <- plot_df %>%
  mutate(
    pathway_clean = factor(pathway_clean, levels = pathway_order),
    mono_y_label = factor(
      mono_y_label,
      levels = rev(annotation_df$mono_y_label)
    )
  )

plot_input_file <- file.path(
  "files",
  "GSEA_Hallmark_Mono_1_to_Mono_6_BubblePlot_Input.csv"
)

write.csv(plot_df, plot_input_file, row.names = FALSE)


# ------------------------------------------------------------------------------
# 7. 绘制 Hallmark GSEA 气泡图
# ------------------------------------------------------------------------------

print("🚀 步骤6: 绘制 Hallmark GSEA 气泡图...")

p_gsea_bubble <- ggplot(
  plot_df,
  aes(x = pathway_clean, y = mono_y_label)
) +
  geom_point(
    aes(size = plot_size, fill = NES),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.90
  ) +
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey85"
  ) +
  scale_size_continuous(
    name = "-log10(FDR)",
    range = c(0.8, 7)
  ) +
  labs(
    title = "Hallmark GSEA of Monocyte Subclusters",
    x = "Hallmark pathway",
    y = "Mono cluster annotation"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 10
    ),
    axis.text.y = element_text(size = 11),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

png_file <- file.path(
  "pictures",
  "GSEA_Hallmark_Mono_1_to_Mono_6_BubblePlot.png"
)

pdf_file <- file.path(
  "pictures",
  "GSEA_Hallmark_Mono_1_to_Mono_6_BubblePlot.pdf"
)

ggsave(
  filename = png_file,
  plot = p_gsea_bubble,
  width = 22,
  height = 6.5,
  dpi = 300,
  limitsize = FALSE
)

ggsave(
  filename = pdf_file,
  plot = p_gsea_bubble,
  width = 22,
  height = 6.5,
  limitsize = FALSE
)

print(paste0("  ✅ 气泡图 PNG 已保存至: ", png_file))
print(paste0("  ✅ 气泡图 PDF 已保存至: ", pdf_file))


# ------------------------------------------------------------------------------
# 8. 额外导出每个 Mono cluster 的 Top enriched / depleted Hallmark 通路
# ------------------------------------------------------------------------------

print("🚀 步骤7: 导出每个 Mono cluster 的 Top Hallmark 通路...")

top_pathways <- gsea_combined %>%
  group_by(mono_cluster_id, mono_annotation) %>%
  arrange(p.adjust, desc(abs(NES)), .by_group = TRUE) %>%
  slice_head(n = 10) %>%
  ungroup()

top_file <- file.path(
  "files",
  "GSEA_Hallmark_Mono_1_to_Mono_6_Top10_PerCluster.csv"
)

write.csv(top_pathways, top_file, row.names = FALSE)

print(paste0("  ✅ Top10 Hallmark 通路已导出至: ", top_file))

print("✅ Script 5 执行完毕：Hallmark GSEA 分析和气泡图均已完成！")