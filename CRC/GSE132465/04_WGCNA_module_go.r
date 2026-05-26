# ==============================================================================
# 单细胞 WGCNA 下游分析: 多模块与整体 GO 富集分析及基因导出 (Signed)
# ==============================================================================

library(clusterProfiler)
library(org.Hs.eg.db)   # 人类注释包；小鼠改成 org.Mm.eg.db
library(ggplot2)
library(qs)

Sys.setenv(http_proxy  = "http://127.0.0.1:7890")
Sys.setenv(https_proxy = "http://127.0.0.1:7890")

WORK_DIR <- '/mnt/disk1/qiuzerui/downloads/CRC/GSE132465'
FILES_DIR <- "files"
PLOTS_DIR <- "plots"
QS_DIR <- "qs"
TARGET_MODULES_DIR <- file.path(FILES_DIR, "Target_Modules_Genes")

setwd(WORK_DIR)

if (!dir.exists(FILES_DIR)) dir.create(FILES_DIR)
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)
if (!dir.exists(TARGET_MODULES_DIR)) dir.create(TARGET_MODULES_DIR)

cat("正在读取 WGCNA 分析结果...\n")

wgcna_results <- qread(file.path(QS_DIR, 'WGCNA_Final_Results_signed.qs'))
gene_module_df <- wgcna_results$gene_module_df

TARGET_COLORS <- c("brown", "blue", "turquoise", "red", "green")
target_gene_df <- gene_module_df[gene_module_df$Module %in% TARGET_COLORS, ]

# ==============================================================================
# 0. 导出目标模块基因
# ==============================================================================
cat("\n正在提取特定目标模块的基因并保存为 CSV...\n")

write.csv(
  target_gene_df,
  file = file.path(FILES_DIR, "Target_Selected_Modules_Genes_signed.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (color in TARGET_COLORS) {
  single_mod_genes <- gene_module_df$Gene[gene_module_df$Module == color]
  
  if (length(single_mod_genes) > 0) {
    output_path <- file.path(TARGET_MODULES_DIR, paste0("Module_", color, "_signed_genes.csv"))
    output_df <- data.frame(Gene = single_mod_genes)
    write.csv(output_df, file = output_path, row.names = FALSE, quote = FALSE)
  } else {
    warning(paste("警告：当前分析结果中未检测到", color, "模块的基因！"))
  }
}

cat("目标模块基因导出完成。\n")

# ==============================================================================
# 1. 基因 ID 转换 Symbol -> Entrez ID
# ==============================================================================
cat("\n开始进行 GO 富集分析的基因 ID 转换...\n")

gene_ids <- bitr(
  target_gene_df$Gene,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

m_gene_df <- merge(target_gene_df, gene_ids, by.x = "Gene", by.y = "SYMBOL")

# ==============================================================================
# 2. 多模块 GO 富集分析 compareCluster
# ==============================================================================
cat("\n正在运行多模块 GO 富集 compareCluster...\n")

go_compare <- compareCluster(
  ENTREZID ~ Module,
  data = m_gene_df,
  fun = "enrichGO",
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",             # BP 生物过程；可改为 MF / CC / ALL
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_compare_df <- as.data.frame(go_compare)

if (nrow(go_compare_df) > 0) {
  
  write.csv(
    go_compare_df,
    file = file.path(FILES_DIR, "Target_Modules_compareCluster_GO_BP_signed.csv"),
    row.names = FALSE
  )
  
  cat("多模块 GO 富集汇总结果已保存。\n")
  
  for (mod in unique(go_compare_df$Cluster)) {
    mod_df <- go_compare_df[go_compare_df$Cluster == mod, ]
    output_name <- paste0("Module_", mod, "_GO_BP_Results_signed.csv")
    
    write.csv(
      mod_df,
      file = file.path(FILES_DIR, output_name),
      row.names = FALSE
    )
  }
  
  cat("各模块独立 GO 富集结果已保存。\n")
  
  p <- dotplot(
    go_compare,
    showCategory = 5,
    title = "GO Biological Process Enrichment Across Modules (Signed)"
  ) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(PLOTS_DIR, "Target_Modules_GO_BP_Dotplot_signed.png"),
    plot = p,
    width = 9,
    height = 11,
    dpi = 300
  )
  
  cat("GO 多模块气泡图已保存。\n")
  
  # ==============================================================================
  # 2.4 提取各模块命中 GO 条目的唯一核心基因
  # ==============================================================================
  cat("\n正在提取各模块命中 GO 条目的核心基因集...\n")
  
  hit_genes_list <- list()
  
  for (mod in unique(go_compare_df$Cluster)) {
    genes_in_terms <- go_compare_df$geneID[go_compare_df$Cluster == mod]
    
    if (length(genes_in_terms) > 0) {
      split_genes <- unlist(strsplit(as.character(genes_in_terms), "/"))
      unique_genes <- unique(trimws(split_genes))
      
      hit_genes_list[[mod]] <- data.frame(
        Module = mod,
        GO_Hit_Gene = unique_genes,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if (length(hit_genes_list) > 0) {
    final_hit_genes_df <- do.call(rbind, hit_genes_list)
    
    hit_output_path <- file.path(
      FILES_DIR,
      "Target_Modules_GO_BP_Hit_Unique_Genes_signed.csv"
    )
    
    write.csv(
      final_hit_genes_df,
      file = hit_output_path,
      row.names = FALSE,
      quote = FALSE
    )
    
    cat("各模块命中 GO 条目的唯一基因列表已保存至:", hit_output_path, "\n")
  }
  
} else {
  warning("未发现任何显著 GO 富集结果。")
}

# ==============================================================================
# 3. 所有目标模块合并后的整体 GO 富集分析
# ==============================================================================
cat("\n正在对所有目标模块基因合并后进行整体 GO 富集分析...\n")

all_target_entrez <- unique(m_gene_df$ENTREZID)

go_all <- enrichGO(
  gene = all_target_entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

go_all_df <- as.data.frame(go_all)

if (nrow(go_all_df) > 0) {
  
  write.csv(
    go_all_df,
    file = file.path(FILES_DIR, "Combined_All_Target_Modules_GO_BP_Results_signed.csv"),
    row.names = FALSE
  )
  
  cat("合并基因的整体 GO 富集结果已保存。\n")
  
  p_all <- dotplot(
    go_all,
    showCategory = 15,
    title = "GO Biological Process Enrichment: All Target Modules Combined"
  ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    file.path(PLOTS_DIR, "Combined_All_Target_Modules_GO_BP_Dotplot_signed.png"),
    plot = p_all,
    width = 8,
    height = 9,
    dpi = 300
  )
  
  cat("合并基因的 GO 气泡图已保存。\n")
  
} else {
  warning("合并后的整体基因集未发现显著 GO 富集结果。")
}

write.csv(target_gene_df, "target_gene_df.csv", row.names = FALSE)

cat("\n============================================\n")
cat("全部 GO 富集分析及结果导出流程运行完毕！\n")
cat("============================================\n")