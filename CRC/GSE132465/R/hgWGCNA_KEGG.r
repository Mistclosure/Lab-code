# ==============================================================================
# 单细胞 WGCNA 下游分析：整体 KEGG + 各模块分别 KEGG
# ==============================================================================

# 如果未安装，请先取消注释并运行：
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)   # 人类注释包；小鼠请改为 org.Mm.eg.db
library(enrichplot)
library(ggplot2)
library(qs)

# ==============================================================================
# 0. 目录配置与数据读取 --------------------------------------------------------
# ==============================================================================

# 如需代理，可取消注释
Sys.setenv(http_proxy  = "http://127.0.0.1:7890")
Sys.setenv(https_proxy = "http://127.0.0.1:7890")

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"
FILES_DIR <- "files"
PLOTS_DIR <- "plots"
QS_DIR <- "qs"

setwd(WORK_DIR)

if (!dir.exists(FILES_DIR)) dir.create(FILES_DIR)
if (!dir.exists(PLOTS_DIR)) dir.create(PLOTS_DIR)

cat("正在读取 WGCNA 分析结果...\n")

wgcna_results <- qread(file.path(QS_DIR, "hdWGCNA_Final_Results_signed.qs"))
gene_module_df <- wgcna_results$gene_module_df

# ==============================================================================
# 1. 设置 KEGG 分析的模块范围 ---------------------------------------------------
# ==============================================================================

# 只分析这些目标模块
TARGET_COLORS <- NULL

# 如果你想分析所有非 grey 模块，可以改成：
# TARGET_COLORS <- NULL

if (is.null(TARGET_COLORS)) {
  kegg_input_df <- gene_module_df[gene_module_df$Module != "grey", ]
  analysis_label <- "All_NonGrey_Modules"
} else {
  kegg_input_df <- gene_module_df[gene_module_df$Module %in% TARGET_COLORS, ]
  analysis_label <- "Target_Modules"
}

if (nrow(kegg_input_df) == 0) {
  stop("没有找到用于 KEGG 分析的基因，请检查 TARGET_COLORS 或 gene_module_df$Module。")
}

cat("用于 KEGG 分析的基因数：", length(unique(kegg_input_df$Gene)), "\n")
cat("涉及模块：", paste(unique(kegg_input_df$Module), collapse = ", "), "\n")

# ==============================================================================
# 2. 基因 ID 转换：Symbol -> Entrez ID ------------------------------------------
# ==============================================================================

cat("\n开始进行 KEGG 富集分析的基因 ID 转换...\n")

gene_ids <- bitr(
  unique(kegg_input_df$Gene),
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

kegg_gene_df <- merge(
  kegg_input_df,
  gene_ids,
  by.x = "Gene",
  by.y = "SYMBOL"
)

kegg_gene_df <- kegg_gene_df[!is.na(kegg_gene_df$ENTREZID), ]

cat("成功转换为 Entrez ID 的基因数：", length(unique(kegg_gene_df$ENTREZID)), "\n")

if (nrow(kegg_gene_df) == 0) {
  stop("没有任何基因成功转换为 Entrez ID，无法进行 KEGG 分析。")
}

# 可选：保存 ID 转换表，便于检查
write.csv(
  kegg_gene_df,
  file = file.path(FILES_DIR, paste0(analysis_label, "_KEGG_Input_Gene_ENTREZID_Map_signed.csv")),
  row.names = FALSE,
  quote = FALSE
)

# 判断 enrichResult 是否有结果
has_kegg_result <- function(x) {
  !is.null(x) && nrow(as.data.frame(x)) > 0
}

# 文件名安全处理
safe_filename <- function(x) {
  gsub("[^A-Za-z0-9_\\-]+", "_", x)
}

# ==============================================================================
# 3. 整体 KEGG：所有目标模块基因合并 -------------------------------------------
# ==============================================================================

cat("\n正在进行整体 KEGG 富集分析：所有目标模块基因合并...\n")

all_entrez <- unique(kegg_gene_df$ENTREZID)

kegg_all <- enrichKEGG(
  gene = all_entrez,
  organism = "hsa",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"
)

if (has_kegg_result(kegg_all)) {
  
  kegg_all <- setReadable(
    kegg_all,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID"
  )
  
  kegg_all_df <- as.data.frame(kegg_all)
  
  write.csv(
    kegg_all_df,
    file = file.path(FILES_DIR, paste0(analysis_label, "_Combined_KEGG_Results_signed.csv")),
    row.names = FALSE
  )
  
  cat("整体 KEGG 结果已保存至 files/ 文件夹。\n")
  
  p_all <- dotplot(
    kegg_all,
    showCategory = 15,
    title = "KEGG Enrichment: Combined Genes"
  ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  ggsave(
    filename = file.path(PLOTS_DIR, paste0(analysis_label, "_Combined_KEGG_Dotplot_signed.png")),
    plot = p_all,
    width = 8,
    height = 9,
    dpi = 300
  )
  
  cat("整体 KEGG 气泡图已保存至 plots/ 文件夹。\n")
  
} else {
  warning("整体基因集未发现显著富集的 KEGG 通路。")
}

# ==============================================================================
# 4. 分别 KEGG：每个模块单独做 KEGG --------------------------------------------
# ==============================================================================

cat("\n正在进行各模块分别 KEGG 富集分析...\n")

modules_to_run <- unique(kegg_gene_df$Module)

separate_kegg_list <- list()

for (mod in modules_to_run) {
  
  cat("\n正在分析模块：", mod, "\n")
  
  mod_entrez <- unique(kegg_gene_df$ENTREZID[kegg_gene_df$Module == mod])
  
  if (length(mod_entrez) < 5) {
    warning(paste0("模块 ", mod, " 可用于 KEGG 的基因数少于 5，跳过该模块。"))
    next
  }
  
  kegg_mod <- enrichKEGG(
    gene = mod_entrez,
    organism = "hsa",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    pAdjustMethod = "BH"
  )
  
  if (has_kegg_result(kegg_mod)) {
    
    kegg_mod <- setReadable(
      kegg_mod,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID"
    )
    
    kegg_mod_df <- as.data.frame(kegg_mod)
    kegg_mod_df$Module <- mod
    
    separate_kegg_list[[mod]] <- kegg_mod_df
    
    mod_safe <- safe_filename(mod)
    
    write.csv(
      kegg_mod_df,
      file = file.path(FILES_DIR, paste0("Module_", mod_safe, "_KEGG_Results_signed.csv")),
      row.names = FALSE
    )
    
    p_mod <- dotplot(
      kegg_mod,
      showCategory = 15,
      title = paste0("KEGG Enrichment: Module ", mod)
    ) +
      theme_bw() +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold")
      )
    
    ggsave(
      filename = file.path(PLOTS_DIR, paste0("Module_", mod_safe, "_KEGG_Dotplot_signed.png")),
      plot = p_mod,
      width = 8,
      height = 9,
      dpi = 300
    )
    
    cat("模块 ", mod, " 的 KEGG 结果和气泡图已保存。\n")
    
  } else {
    warning(paste0("模块 ", mod, " 未发现显著富集的 KEGG 通路。"))
  }
}

# 合并保存所有模块分别 KEGG 的结果
if (length(separate_kegg_list) > 0) {
  
  separate_kegg_df <- do.call(rbind, separate_kegg_list)
  
  write.csv(
    separate_kegg_df,
    file = file.path(FILES_DIR, paste0(analysis_label, "_Separate_Modules_KEGG_All_Results_signed.csv")),
    row.names = FALSE
  )
  
  cat("\n所有模块分别 KEGG 的合并结果已保存至 files/ 文件夹。\n")
  
} else {
  warning("没有任何模块得到显著 KEGG 富集结果。")
}

cat("\n============================================\n")
cat("整体 KEGG 和各模块分别 KEGG 分析完成！\n")
cat("============================================\n")