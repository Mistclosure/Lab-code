# ==============================================================================
# 单细胞 WGCNA 下游分析: 多模块与整体 KEGG 富集分析及基因导出 (Signed)
# ==============================================================================

# 如果未安装，请先取消注释并运行以下 BiocManager 安装命令：
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "ggplot2"))

library(clusterProfiler)
library(org.Hs.eg.db) # 人类注释包。如果是小鼠，请使用 org.Mm.eg.db
library(ggplot2)
library(qs)

# ==============================================================================
# 0. 目录配置与数据读取 --------------------------------------------------------
# ==============================================================================
# 设置为vpn端口
Sys.setenv(http_proxy  = "http://127.0.0.1:7890")
Sys.setenv(https_proxy = "http://127.0.0.1:7890")

WORK_DIR <- '/mnt/disk1/qiuzerui/downloads/CRC/GSE132465'
FILES_DIR <- "files"
PLOTS_DIR <- "plots"
QS_DIR <- "qs"
TARGET_MODULES_DIR <- file.path(FILES_DIR, "Target_Modules_Genes") # 新增：存放按模块独立导出的基因文件目录

setwd(WORK_DIR)

# 确保输出文件夹存在
if (!dir.exists(FILES_DIR)) dir.create(FILES_DIR)
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)

cat("正在读取 WGCNA 分析结果...\n")
# 读取之前保存的 qs 结果列表，并提取基因模块映射表
wgcna_results <- qread(file.path(QS_DIR, 'WGCNA_Final_Results_signed.qs'))
gene_module_df <- wgcna_results$gene_module_df

# 提取你关心的目标模块
TARGET_COLORS <- c("brown", "blue", "turquoise", "red", "green")
target_gene_df <- gene_module_df[gene_module_df$Module %in% TARGET_COLORS, ]

# ==============================================================================
# [新增模块] 导出：特定目标模块的基因到 CSV 文件
# ==============================================================================
cat("\n正在提取特定目标模块的基因并保存为CSV...\n")

# 1. 过滤核心模块并保存总文件
write.csv(target_gene_df, file = file.path(FILES_DIR, "Target_Selected_Modules_Genes_signed.csv"), row.names = FALSE, quote = FALSE)
cat("已将目标模块的基因合并保存至 ", file.path(FILES_DIR, "Target_Selected_Modules_Genes_signed.csv"), "\n")

# 2. 循环保存单模块文件 (输出为 CSV)
if (!dir.exists(TARGET_MODULES_DIR)) dir.create(TARGET_MODULES_DIR)

for (color in TARGET_COLORS) {
  single_mod_genes <- gene_module_df$Gene[gene_module_df$Module == color]
  if (length(single_mod_genes) > 0) {
    output_path <- file.path(TARGET_MODULES_DIR, paste0("Module_", color, "_signed_genes.csv"))
    # 转换为数据框保证带列名
    output_df <- data.frame(Gene = single_mod_genes)
    write.csv(output_df, file = output_path, row.names = FALSE, quote = FALSE)
  } else {
    warning(paste("警告：当前分析结果中未检测到", color, "模块的基因！"))
  }
}
cat("已将各目标模块的基因列表分别独立保存至目录：", TARGET_MODULES_DIR, "/ 下\n\n")

# ==============================================================================
# 1. 基因 ID 转换 (Symbol -> Entrez ID) ----------------------------------------
# ==============================================================================
cat("开始进行 KEGG 富集分析的基因 ID 转换...\n")
# KEGG 数据库底层主要识别 Entrez ID，因此需要先进行转换
gene_ids <- bitr(target_gene_df$Gene,
                 fromType = "SYMBOL",   # 原本的基因名类型
                 toType = "ENTREZID",
                 OrgDb = org.Hs.eg.db)  # 人类数据库

# 将转换后的 Entrez ID 合并回我们原有的模块分类数据框中
m_gene_df <- merge(target_gene_df, gene_ids, by.x = "Gene", by.y = "SYMBOL")

# ==============================================================================
# 2. 多模块整体 KEGG 富集分析 (compareCluster) ---------------------------------
# ==============================================================================
cat("\n正在运行多模块比对 KEGG 富集 (compareCluster)...\n")
# 公式 ENTREZID ~ Module 意为：按 Module 分组，对 ENTREZID 进行富集
kegg_compare <- compareCluster(ENTREZID ~ Module, 
                               data = m_gene_df, 
                               fun = "enrichKEGG", 
                               organism = "hsa",    # 人类是 "hsa"
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.2)

if (!is.null(kegg_compare)) {
  # 将结果中的 Entrez ID 转换回可读的 Gene Symbol
  kegg_compare <- setReadable(kegg_compare, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # ----------------------------------------------------------------------------
  # 2.1 保存完整的多模块富集表格数据
  # ----------------------------------------------------------------------------
  kegg_compare_df <- as.data.frame(kegg_compare)
  write.csv(kegg_compare_df, 
            file = file.path(FILES_DIR, "Target_Modules_compareCluster_KEGG_signed.csv"), 
            row.names = FALSE)
  cat("所有模块 compareCluster 汇总结果已保存至 files/ 文件夹。\n")
  
  # ----------------------------------------------------------------------------
  # 2.2 循环输出每一个模块单独的 KEGG 结果
  # ----------------------------------------------------------------------------
  cat("正在提取并保存各个独立模块的 KEGG 结果...\n")
  for (mod in unique(kegg_compare_df$Cluster)) {
    mod_df <- kegg_compare_df[kegg_compare_df$Cluster == mod, ]
    output_name <- paste0("Module_", mod, "_KEGG_Results_signed.csv")
    write.csv(mod_df, 
              file = file.path(FILES_DIR, output_name), 
              row.names = FALSE)
  }
  cat("各模块独立的 KEGG 结果已拆分并保存至 files/ 文件夹。\n")
  
  # ----------------------------------------------------------------------------
  # 2.3 多模块横向对比可视化（高级气泡图）
  # ----------------------------------------------------------------------------
  # showCategory = 5 代表每个模块展示前 5 个最显著的通路
  p <- dotplot(kegg_compare, showCategory = 5, title = "KEGG Pathway Enrichment Across Modules (Signed)") +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
      axis.text.y = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 保存图片到 plots 文件夹
  ggsave(file.path(PLOTS_DIR, "Target_Modules_KEGG_Dotplot_signed.png"), plot = p, width = 9, height = 11, dpi = 300)
  cat("KEGG 对比气泡图已成功绘制并保存至 plots/ 文件夹！\n")
  
  # ----------------------------------------------------------------------------
  # 2.4 [新增功能] 提取并保存各模块命中 KEGG 的唯一核心基因列表
  # ----------------------------------------------------------------------------
  cat("\n正在提取各模块命中 KEGG 通路的核心基因集...\n")
  hit_genes_list <- list()
  
  for (mod in unique(kegg_compare_df$Cluster)) {
    # 提取该模块在所有显著富集通路中的 geneID 列
    genes_in_pathways <- kegg_compare_df$geneID[kegg_compare_df$Cluster == mod]
    
    if (length(genes_in_pathways) > 0) {
      # 按照 "/" 拆分字符串，将所有的基因打散成向量
      split_genes <- unlist(strsplit(as.character(genes_in_pathways), "/"))
      # 取唯一值去重
      unique_genes <- unique(split_genes)
      unique_genes <- trimws(unique_genes)
      
      # 保存为数据框格式
      hit_genes_list[[mod]] <- data.frame(
        Module = mod,
        KEGG_Hit_Gene = unique_genes,
        stringsAsFactors = FALSE
      )
    }
  }
  
  # 将所有模块的结果合并为一个总的数据框
  if (length(hit_genes_list) > 0) {
    final_hit_genes_df <- do.call(rbind, hit_genes_list)
    # 导出为 CSV
    hit_output_path <- file.path(FILES_DIR, "Target_Modules_KEGG_Hit_Unique_Genes_signed.csv")
    write.csv(final_hit_genes_df, file = hit_output_path, row.names = FALSE, quote = FALSE)
    cat(">>> 各模块命中 KEGG 通路的唯一基因列表已成功保存至:", hit_output_path, "<<<\n")
  } else {
    cat("提示：未提取到任何命中 KEGG 的基因。\n")
  }

} else {
  warning("未发现任何在设定阈值下显著富集的多模块 KEGG 通路。")
}

# ==============================================================================
# 3. 所有目标模块基因合并的整体 KEGG 富集分析 ----------------------------------
# ==============================================================================
cat("\n正在对所有目标模块基因合并后进行整体 KEGG 富集分析 (enrichKEGG)...\n")

# 提取所有目标模块去重后的 Entrez ID
all_target_entrez <- unique(m_gene_df$ENTREZID)

kegg_all <- enrichKEGG(gene = all_target_entrez, 
                       organism = "hsa", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.2)

if (!is.null(kegg_all)) {
  # 转换回可读的 Gene Symbol
  kegg_all <- setReadable(kegg_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  # 保存合并后的整体 KEGG 富集结果
  write.csv(as.data.frame(kegg_all), 
            file = file.path(FILES_DIR, "Combined_All_Target_Modules_KEGG_Results_signed.csv"), 
            row.names = FALSE)
  cat("合并基因的整体 KEGG 富集结果已保存至 files/ 文件夹。\n")
  
  # 绘制合并后的整体气泡图
  p_all <- dotplot(kegg_all, showCategory = 15, title = "KEGG Enrichment (All Target Modules Combined - Signed)") +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(file.path(PLOTS_DIR, "Combined_All_Target_Modules_KEGG_Dotplot_signed.png"), plot = p_all, width = 8, height = 9, dpi = 300)
  cat("合并基因的 KEGG 气泡图已保存至 plots/ 文件夹！\n")
  
} else {
  warning("合并后的整体基因集未发现显著富集的 KEGG 通路。")
}

cat("\n============================================\n")
cat("全部 KEGG 分析及结果导出流程运行完毕！\n")
cat("============================================\n")

write.csv(target_gene_df, "target_gene_df.csv", row.names = FALSE)