# ==============================================================================
# cNMF program Top200 genes KEGG 富集分析
# 基于 CiliaHub 限定基因集的 cNMF spectra / marker ranking 结果
# ==============================================================================

# 如果未安装，请先取消注释并运行：
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)   # 人类注释包
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# ==============================================================================
# 0. 目录配置与参数 ------------------------------------------------------------
# ==============================================================================

# 如需代理，可取消注释
# Sys.setenv(http_proxy  = "http://127.0.0.1:7890")
# Sys.setenv(https_proxy = "http://127.0.0.1:7890")

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"
setwd(WORK_DIR)

# 输入文件：
# 1) 默认优先使用 marker ranking long format:
#    program, rank, gene, spectra_score
# 2) 如果没有该文件，也支持 Program x Gene matrix:
#    consensus_programs_by_gene.k_7.csv
CNMF_TABLE_PATH <- file.path(
  WORK_DIR,
  "cNMF_CiliaCarta_GSE132465",
  "tables",
  "marker_gene_rankings_by_program.k_5.csv"
)

CNMF_MATRIX_FALLBACK_PATH <- file.path(
  WORK_DIR,
  "cNMF_CiliaCarta_GSE132465",
  "tables",
  "consensus_programs_by_gene.k_5.csv"
)

OUT_DIR <- file.path(WORK_DIR, "KEGG_cillacorta_cNMF_program_top200")
FILES_DIR <- file.path(OUT_DIR, "files")
PLOTS_DIR <- file.path(OUT_DIR, "plots")
PROGRAM_GENES_DIR <- file.path(FILES_DIR, "Program_Top200_Genes")
PROGRAM_KEGG_DIR <- file.path(FILES_DIR, "Program_KEGG_Results")

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FILES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PROGRAM_GENES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PROGRAM_KEGG_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N <- 200
MIN_ENTREZ_GENES <- 5
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.2
PADJUST_METHOD <- "BH"
SHOW_CATEGORY_PER_PROGRAM <- 15
SHOW_CATEGORY_COMBINED_PER_PROGRAM <- 8

cat("WORK_DIR:", WORK_DIR, "\n")
cat("CNMF_TABLE_PATH:", CNMF_TABLE_PATH, "\n")
cat("CNMF_MATRIX_FALLBACK_PATH:", CNMF_MATRIX_FALLBACK_PATH, "\n")
cat("OUT_DIR:", OUT_DIR, "\n")

# ==============================================================================
# 1. 工具函数 ------------------------------------------------------------------
# ==============================================================================

has_kegg_result <- function(x) {
  !is.null(x) && nrow(as.data.frame(x)) > 0
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

gene_ratio_to_numeric <- function(x) {
  vapply(
    strsplit(as.character(x), "/"),
    function(z) {
      if (length(z) != 2) return(NA_real_)
      as.numeric(z[1]) / as.numeric(z[2])
    },
    numeric(1)
  )
}

write_empty_kegg_file <- function(program, reason) {
  empty_df <- data.frame(
    program = program,
    message = reason,
    stringsAsFactors = FALSE
  )
  out_file <- file.path(
    PROGRAM_KEGG_DIR,
    paste0(safe_filename(program), "_KEGG_Results_top200.csv")
  )
  write.csv(empty_df, file = out_file, row.names = FALSE, quote = FALSE)
  
  msg_file <- file.path(
    PROGRAM_KEGG_DIR,
    paste0(safe_filename(program), "_KEGG_no_significant_result.txt")
  )
  writeLines(reason, con = msg_file)
}

save_plot_pdf_png <- function(plot_obj, file_prefix, width, height) {
  ggsave(
    filename = paste0(file_prefix, ".png"),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    filename = paste0(file_prefix, ".pdf"),
    plot = plot_obj,
    width = width,
    height = height
  )
}

# ==============================================================================
# 2. 读取 cNMF 结果并提取每个 Program Top200 genes -----------------------------
# ==============================================================================

cat("\n正在读取 cNMF program gene score 表...\n")

if (file.exists(CNMF_TABLE_PATH)) {
  
  cat("检测到 marker ranking long format，使用该文件：", CNMF_TABLE_PATH, "\n")
  
  cnmf_long <- read.csv(
    CNMF_TABLE_PATH,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  required_cols <- c("program", "gene", "spectra_score")
  if (!all(required_cols %in% colnames(cnmf_long))) {
    stop(
      "marker ranking 表必须包含列：",
      paste(required_cols, collapse = ", "),
      "\n当前列名：",
      paste(colnames(cnmf_long), collapse = ", ")
    )
  }
  
  cnmf_long <- cnmf_long %>%
    mutate(
      program = as.character(program),
      gene = as.character(gene),
      spectra_score = as.numeric(spectra_score)
    ) %>%
    filter(!is.na(program), program != "", !is.na(gene), gene != "", !is.na(spectra_score))
  
} else if (file.exists(CNMF_MATRIX_FALLBACK_PATH)) {
  
  cat("未找到 marker ranking 表，改用 Program x Gene matrix：", CNMF_MATRIX_FALLBACK_PATH, "\n")
  
  cnmf_mat <- read.csv(
    CNMF_MATRIX_FALLBACK_PATH,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  cnmf_mat$program <- rownames(cnmf_mat)
  
  cnmf_long <- cnmf_mat %>%
    pivot_longer(
      cols = -program,
      names_to = "gene",
      values_to = "spectra_score"
    ) %>%
    mutate(
      program = as.character(program),
      gene = as.character(gene),
      spectra_score = as.numeric(spectra_score)
    ) %>%
    filter(!is.na(program), program != "", !is.na(gene), gene != "", !is.na(spectra_score))
  
} else {
  stop(
    "找不到 cNMF 输入表。\n",
    "请检查 CNMF_TABLE_PATH 或 CNMF_MATRIX_FALLBACK_PATH。\n",
    "当前 CNMF_TABLE_PATH: ", CNMF_TABLE_PATH, "\n",
    "当前 CNMF_MATRIX_FALLBACK_PATH: ", CNMF_MATRIX_FALLBACK_PATH
  )
}

program_order <- unique(cnmf_long$program)
program_order <- program_order[order(as.integer(gsub("[^0-9]+", "", program_order)))]

cat("检测到 Program：", paste(program_order, collapse = ", "), "\n")

top200_df <- cnmf_long %>%
  group_by(program) %>%
  arrange(desc(spectra_score), .by_group = TRUE) %>%
  slice_head(n = TOP_N) %>%
  mutate(top_rank = row_number()) %>%
  ungroup() %>%
  select(program, top_rank, gene, spectra_score)

top200_all_file <- file.path(FILES_DIR, "program_top200_genes_all.csv")
write.csv(top200_df, file = top200_all_file, row.names = FALSE, quote = FALSE)
cat("Top200 genes 总表已保存：", top200_all_file, "\n")

for (program in program_order) {
  program_df <- top200_df[top200_df$program == program, ]
  program_file <- file.path(
    PROGRAM_GENES_DIR,
    paste0(safe_filename(program), "_top200_genes.csv")
  )
  write.csv(program_df, file = program_file, row.names = FALSE, quote = FALSE)
}
cat("每个 Program 的 Top200 genes 已分别保存至：", PROGRAM_GENES_DIR, "\n")

# ==============================================================================
# 3. 对每个 Program 分别进行 KEGG 富集 -----------------------------------------
# ==============================================================================

cat("\n开始逐个 Program 进行 KEGG 富集分析...\n")

program_summary_list <- list()
kegg_all_list <- list()
id_map_list <- list()
unmapped_list <- list()

for (program in program_order) {
  
  cat("\n正在分析：", program, "\n")
  
  program_top_df <- top200_df[top200_df$program == program, ]
  input_genes <- unique(program_top_df$gene)
  input_gene_count <- length(input_genes)
  
  if (input_gene_count == 0) {
    warning(paste0(program, " 没有可用于 KEGG 的 Top genes，跳过。"))
    write_empty_kegg_file(program, "No input genes.")
    next
  }
  
  gene_ids <- suppressWarnings(
    bitr(
      input_genes,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db
    )
  )
  
  if (is.null(gene_ids) || nrow(gene_ids) == 0) {
    gene_ids <- data.frame(SYMBOL = character(), ENTREZID = character())
  }
  
  gene_ids <- gene_ids[!duplicated(gene_ids[, c("SYMBOL", "ENTREZID")]), ]
  mapped_symbols <- unique(gene_ids$SYMBOL)
  unmapped_symbols <- setdiff(input_genes, mapped_symbols)
  entrez_ids <- unique(gene_ids$ENTREZID)
  mapped_count <- length(unique(mapped_symbols))
  
  if (nrow(gene_ids) > 0) {
    id_map_df <- data.frame(
      program = rep(program, nrow(gene_ids)),
      SYMBOL = gene_ids$SYMBOL,
      ENTREZID = gene_ids$ENTREZID,
      stringsAsFactors = FALSE
    )
  } else {
    id_map_df <- data.frame(
      program = character(),
      SYMBOL = character(),
      ENTREZID = character(),
      stringsAsFactors = FALSE
    )
  }
  id_map_list[[program]] <- id_map_df
  
  if (length(unmapped_symbols) > 0) {
    unmapped_df <- data.frame(
      program = rep(program, length(unmapped_symbols)),
      unmapped_gene_symbol = unmapped_symbols,
      stringsAsFactors = FALSE
    )
  } else {
    unmapped_df <- data.frame(
      program = character(),
      unmapped_gene_symbol = character(),
      stringsAsFactors = FALSE
    )
  }
  unmapped_list[[program]] <- unmapped_df
  
  cat("Top genes:", input_gene_count, "\n")
  cat("Entrez ID 转换成功 gene symbol 数:", mapped_count, "\n")
  cat("未转换 gene symbol 数:", length(unmapped_symbols), "\n")
  
  if (length(entrez_ids) < MIN_ENTREZ_GENES) {
    reason <- paste0(
      program,
      " 转换到 Entrez ID 的基因数少于 ",
      MIN_ENTREZ_GENES,
      "，跳过 KEGG。"
    )
    warning(reason)
    write_empty_kegg_file(program, reason)
    
    program_summary_list[[program]] <- data.frame(
      program = program,
      input_gene_count = input_gene_count,
      mapped_gene_symbol_count = mapped_count,
      entrez_id_count = length(entrez_ids),
      unmapped_gene_symbols = paste(unmapped_symbols, collapse = "/"),
      significant_kegg_count = 0,
      stringsAsFactors = FALSE
    )
    next
  }
  
  kegg_res <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = PVALUE_CUTOFF,
    qvalueCutoff = QVALUE_CUTOFF,
    pAdjustMethod = PADJUST_METHOD
  )
  
  if (has_kegg_result(kegg_res)) {
    
    kegg_res <- setReadable(
      kegg_res,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    
    kegg_df <- as.data.frame(kegg_res)
    kegg_df$program <- program
    kegg_df$input_gene_count <- input_gene_count
    kegg_df$mapped_gene_symbol_count <- mapped_count
    kegg_df$entrez_id_count <- length(entrez_ids)
    kegg_df$unmapped_gene_symbols <- paste(unmapped_symbols, collapse = "/")
    kegg_df$pathway_hit_gene_symbols <- kegg_df$geneID
    
    kegg_df <- kegg_df[, c(
      "program",
      "input_gene_count",
      "mapped_gene_symbol_count",
      "entrez_id_count",
      "unmapped_gene_symbols",
      "ID",
      "Description",
      "GeneRatio",
      "BgRatio",
      "pvalue",
      "p.adjust",
      "qvalue",
      "geneID",
      "pathway_hit_gene_symbols",
      "Count"
    )]
    
    kegg_all_list[[program]] <- kegg_df
    
    program_kegg_file <- file.path(
      PROGRAM_KEGG_DIR,
      paste0(safe_filename(program), "_KEGG_Results_top200.csv")
    )
    write.csv(kegg_df, file = program_kegg_file, row.names = FALSE, quote = FALSE)
    
    p_mod <- dotplot(
      kegg_res,
      showCategory = SHOW_CATEGORY_PER_PROGRAM,
      title = paste0("KEGG Enrichment: ", program, " Top200 genes")
    ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    plot_prefix <- file.path(
      PLOTS_DIR,
      paste0(safe_filename(program), "_KEGG_Dotplot_top200")
    )
    save_plot_pdf_png(p_mod, plot_prefix, width = 8, height = 9)
    
    cat(program, " KEGG 结果和气泡图已保存。\n")
    
    significant_count <- nrow(kegg_df)
    
  } else {
    
    reason <- paste0(program, " 未发现显著富集的 KEGG 通路。")
    warning(reason)
    write_empty_kegg_file(program, reason)
    significant_count <- 0
  }
  
  program_summary_list[[program]] <- data.frame(
    program = program,
    input_gene_count = input_gene_count,
    mapped_gene_symbol_count = mapped_count,
    entrez_id_count = length(entrez_ids),
    unmapped_gene_symbols = paste(unmapped_symbols, collapse = "/"),
    significant_kegg_count = significant_count,
    stringsAsFactors = FALSE
  )
}

# ==============================================================================
# 4. 保存 ID 转换、未转换基因、合并 KEGG 结果 ----------------------------------
# ==============================================================================

id_map_all <- if (length(id_map_list) > 0) {
  do.call(rbind, id_map_list)
} else {
  data.frame(program = character(), SYMBOL = character(), ENTREZID = character())
}

unmapped_all <- if (length(unmapped_list) > 0) {
  do.call(rbind, unmapped_list)
} else {
  data.frame(program = character(), unmapped_gene_symbol = character())
}

program_summary_df <- if (length(program_summary_list) > 0) {
  do.call(rbind, program_summary_list)
} else {
  data.frame(
    program = character(),
    input_gene_count = integer(),
    mapped_gene_symbol_count = integer(),
    entrez_id_count = integer(),
    unmapped_gene_symbols = character(),
    significant_kegg_count = integer()
  )
}

write.csv(
  id_map_all,
  file = file.path(FILES_DIR, "program_top200_SYMBOL_to_ENTREZID_map.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  unmapped_all,
  file = file.path(FILES_DIR, "program_top200_unmapped_gene_symbols.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  program_summary_df,
  file = file.path(FILES_DIR, "program_top200_KEGG_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

if (length(kegg_all_list) > 0) {
  kegg_all_df <- do.call(rbind, kegg_all_list)
} else {
  kegg_all_df <- data.frame(
    program = character(),
    input_gene_count = integer(),
    mapped_gene_symbol_count = integer(),
    entrez_id_count = integer(),
    unmapped_gene_symbols = character(),
    ID = character(),
    Description = character(),
    GeneRatio = character(),
    BgRatio = character(),
    pvalue = numeric(),
    p.adjust = numeric(),
    qvalue = numeric(),
    geneID = character(),
    pathway_hit_gene_symbols = character(),
    Count = integer()
  )
}

combined_kegg_file <- file.path(FILES_DIR, "KEGG_enrichment_all_programs_top200.csv")
write.csv(kegg_all_df, file = combined_kegg_file, row.names = FALSE, quote = FALSE)
cat("\n合并 KEGG 结果已保存：", combined_kegg_file, "\n")

# ==============================================================================
# 5. 合并 KEGG 气泡图 -----------------------------------------------------------
# ==============================================================================

if (nrow(kegg_all_df) > 0) {
  
  plot_df <- kegg_all_df %>%
    group_by(program) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = SHOW_CATEGORY_COMBINED_PER_PROGRAM) %>%
    ungroup()
  
  plot_df$GeneRatio_numeric <- gene_ratio_to_numeric(plot_df$GeneRatio)
  plot_df$minus_log10_padj <- -log10(plot_df$p.adjust)
  plot_df$minus_log10_padj[is.infinite(plot_df$minus_log10_padj)] <- max(
    plot_df$minus_log10_padj[is.finite(plot_df$minus_log10_padj)],
    na.rm = TRUE
  )
  
  # pathway 过长时自动换行，同时限制展示通路集合为各 Program 的 top 条目并集
  pathway_order <- plot_df %>%
    group_by(Description) %>%
    summarise(best_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(best_padj) %>%
    pull(Description)
  
  plot_df$Description_wrapped <- str_wrap(plot_df$Description, width = 45)
  description_level_map <- setNames(str_wrap(pathway_order, width = 45), pathway_order)
  plot_df$Description_wrapped <- factor(
    plot_df$Description_wrapped,
    levels = rev(unique(description_level_map))
  )
  
  plot_df$program <- factor(plot_df$program, levels = program_order)
  
  p_compare <- ggplot(
    plot_df,
    aes(
      x = program,
      y = Description_wrapped,
      size = Count,
      color = minus_log10_padj
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "#2C7BB6",
      high = "#D7191C",
      name = "-log10(adj.P)"
    ) +
    scale_size_continuous(name = "Count") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = "KEGG Pathway Enrichment Across cNMF Programs (Top200 genes)",
      x = "cNMF Program",
      y = "KEGG pathway"
    )
  
  compare_prefix <- file.path(PLOTS_DIR, "All_Programs_KEGG_Bubbleplot_top200")
  save_plot_pdf_png(
    p_compare,
    compare_prefix,
    width = max(9, length(program_order) * 0.8),
    height = max(8, length(unique(plot_df$Description)) * 0.28)
  )
  
  cat("合并 KEGG 气泡图已保存至：", PLOTS_DIR, "\n")
  
} else {
  warning("没有任何 Program 得到显著 KEGG 富集结果，跳过合并气泡图。")
  writeLines(
    "No significant KEGG pathway was found for any cNMF program.",
    con = file.path(PLOTS_DIR, "All_Programs_KEGG_Bubbleplot_top200_NOT_CREATED.txt")
  )
}

# ==============================================================================
# 6. 最终汇总 ------------------------------------------------------------------
# ==============================================================================

cat("\n============================================\n")
cat("cNMF Program Top200 genes KEGG 分析完成！\n")
cat("============================================\n")

for (program in program_order) {
  one_summary <- program_summary_df[program_summary_df$program == program, ]
  if (nrow(one_summary) == 0) next
  
  cat(
    program,
    " | Top genes:",
    one_summary$input_gene_count,
    " | Entrez 转换成功:",
    one_summary$mapped_gene_symbol_count,
    " | 显著 KEGG pathway:",
    one_summary$significant_kegg_count,
    "\n"
  )
}

cat("\n所有输出文件所在目录：", OUT_DIR, "\n")
cat("Top200 gene 文件目录：", PROGRAM_GENES_DIR, "\n")
cat("KEGG 表格目录：", FILES_DIR, "\n")
cat("KEGG 图片目录：", PLOTS_DIR, "\n")
