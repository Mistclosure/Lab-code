# ==============================================================================
# 04_Selected_Modules_KEGG.R
# 对指定模块基因进行 KEGG 富集分析并输出气泡图
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
})

# ==============================================================================
# 1. 参数设置
# ==============================================================================

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"
setwd(WORK_DIR)

FILES_DIR <- "files"

# 这里改成你想分析的模块颜色
SELECT_MODULE_COLORS <- c("brown")

# 步骤 13 保存模块基因的目录
SELECT_MODULES_DIR <- file.path(FILES_DIR, "Selected_Modules_By_Color")

# KEGG 输出目录
KEGG_DIR <- file.path(FILES_DIR, "Selected_Modules_KEGG")
dir.create(KEGG_DIR, showWarnings = FALSE, recursive = TRUE)

# 富集参数
P_VALUE_CUTOFF <- 0.05
Q_VALUE_CUTOFF <- 0.2
SHOW_CATEGORY <- 20

# ==============================================================================
# 2. 单个模块 KEGG 富集函数
# ==============================================================================

run_kegg_for_module <- function(module_color) {
  
  cat("\n==============================\n")
  cat("正在分析模块：", module_color, "\n")
  cat("==============================\n")
  
  module_gene_file <- file.path(
    SELECT_MODULES_DIR,
    paste0("Module_", module_color, "_genes.csv")
  )
  
  if (!file.exists(module_gene_file)) {
    warning("未找到模块基因文件：", module_gene_file)
    return(NULL)
  }
  
  module_genes <- read.csv(
    module_gene_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )$Gene
  
  module_genes <- unique(na.omit(as.character(module_genes)))
  
  cat("模块基因数：", length(module_genes), "\n")
  
  if (length(module_genes) < 5) {
    warning("模块 ", module_color, " 基因数过少，跳过 KEGG。")
    return(NULL)
  }
  
  gene_convert <- bitr(
    module_genes,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )
  
  gene_convert <- gene_convert[!duplicated(gene_convert$ENTREZID), ]
  entrez_ids <- unique(gene_convert$ENTREZID)
  
  cat("成功转换 ENTREZID 基因数：", length(entrez_ids), "\n")
  
  write.csv(
    gene_convert,
    file = file.path(KEGG_DIR, paste0("Gene_SYMBOL_to_ENTREZ_", module_color, ".csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  if (length(entrez_ids) < 5) {
    warning("模块 ", module_color, " 可转换 ENTREZID 的基因过少，跳过 KEGG。")
    return(NULL)
  }
  
  kegg_res <- enrichKEGG(
    gene = entrez_ids,
    organism = "hsa",
    pvalueCutoff = P_VALUE_CUTOFF,
    qvalueCutoff = Q_VALUE_CUTOFF
  )
  
  if (is.null(kegg_res) || nrow(as.data.frame(kegg_res)) == 0) {
    warning("模块 ", module_color, " 未富集到显著 KEGG 通路。")
    return(NULL)
  }
  
  kegg_res <- setReadable(
    kegg_res,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  kegg_df <- as.data.frame(kegg_res)
  
  write.csv(
    kegg_df,
    file = file.path(KEGG_DIR, paste0("KEGG_", module_color, "_result.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  p <- dotplot(
    kegg_res,
    showCategory = SHOW_CATEGORY,
    title = paste0("KEGG enrichment: ", module_color, " module")
  ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10)
    )
  
  ggsave(
    filename = file.path(KEGG_DIR, paste0("KEGG_", module_color, "_bubbleplot.png")),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(KEGG_DIR, paste0("KEGG_", module_color, "_bubbleplot.pdf")),
    plot = p,
    width = 8,
    height = 6
  )
  
  cat("模块 ", module_color, " KEGG 富集完成。\n")
  
  return(kegg_df)
}

# ==============================================================================
# 3. 批量运行所有指定模块
# ==============================================================================

all_kegg_results <- list()

for (module_color in SELECT_MODULE_COLORS) {
  all_kegg_results[[module_color]] <- run_kegg_for_module(module_color)
}

# ==============================================================================
# 4. 合并所有模块 KEGG 结果
# ==============================================================================

combined_kegg <- bind_rows(
  lapply(names(all_kegg_results), function(module_color) {
    df <- all_kegg_results[[module_color]]
    if (is.null(df)) return(NULL)
    df$Module <- module_color
    df
  })
)

if (!is.null(combined_kegg) && nrow(combined_kegg) > 0) {
  write.csv(
    combined_kegg,
    file = file.path(KEGG_DIR, "KEGG_all_selected_modules_combined.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

cat("\n全部 KEGG 分析完成。\n")
cat("结果目录：", KEGG_DIR, "\n")