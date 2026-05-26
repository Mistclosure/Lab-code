# ==============================================================================
# Script: Monocytes 整体 Cold vs RT 差异基因分析、火山图及 GO 富集
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 加载必要的 R 包
library(Seurat)
library(tidyverse)
library(qs)
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠基因注释数据库
library(ggrepel)      # 用于火山图基因标签防重叠

# 创建输出文件夹
if (!dir.exists("files")) dir.create("files")
if (!dir.exists("plots")) dir.create("plots")

# ------------------------------------------------------------------------------
# 1. 读取单核细胞精细聚类对象并设置分组
# ------------------------------------------------------------------------------
print("🚀 步骤1: 正在加载 Monocytes 子集对象...")
mono_cells <- qread("pbmc_monocytes_sub-clustered.qs")

# 过滤掉 TN_30C 组（如果数据中包含的话）
mono_cells <- subset(mono_cells, subset = Group == 'TN_30C', invert = TRUE)

# 【关键点】将当前鉴定的标准切换为 Group，以便比较 Cold 和 RT
Idents(mono_cells) <- "Group"

# 注意：请确保您的 metadata 中包含确切的 "Cold" 和 "RT" 标签
# 如果您的命名是 "Cold_4C" 或 "RT_25C"，请在下方 ident.1 和 ident.2 中对应修改！
group_cold <- "Cold_4C" 
group_rt   <- "RT_25C"

# ------------------------------------------------------------------------------
# 2. 计算 Cold vs RT 的所有差异基因
# ------------------------------------------------------------------------------
print(paste0("🚀 步骤2: 正在计算 ", group_cold, " vs ", group_rt, " 的差异基因..."))
# 注意：为了画出完整的火山图，这里将 logfc.threshold 设为 0（计算所有基因）
deg_results <- FindMarkers(mono_cells, 
                           ident.1 = group_cold, 
                           ident.2 = group_rt, 
                           logfc.threshold = 0, 
                           min.pct = 0.1, 
                           verbose = FALSE)

# 整理数据框，添加基因列和上下调标签
# 设定阈值：|log2FC| > 0.25 且 p_val_adj < 0.05
deg_results$gene <- rownames(deg_results)
fc_threshold <- 0.25
p_val_threshold <- 0.05

deg_results <- deg_results %>%
  mutate(Significance = case_when(
    p_val_adj < p_val_threshold & avg_log2FC > fc_threshold ~ "Up in Cold",
    p_val_adj < p_val_threshold & avg_log2FC < -fc_threshold ~ "Down in Cold",
    TRUE ~ "Not Sig"
  ))

# 保存差异分析全表
#write.csv(deg_results, file = file.path("files", "DEG_Monocytes_Cold_vs_RT.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. 绘制差异基因火山图 (Volcano Plot)
# ------------------------------------------------------------------------------
print("🚀 步骤3: 正在绘制并保存火山图...")

# 提取每组 top 10 显著的基因用于在火山图上添加文本标签
top_genes <- deg_results %>%
  filter(Significance != "Not Sig") %>%
  group_by(Significance) %>%
  top_n(n = 10, wt = abs(avg_log2FC)) # 按 log2FC 绝对值取 top10，也可按 -p_val_adj 取

# 处理一下极端 p 值（如果有 p_val_adj 为 0 的，替换为一个极小值防止 -log10 报错为 Inf）
min_nonzero_pval <- min(deg_results$p_val_adj[deg_results$p_val_adj > 0])
deg_results$p_val_adj[deg_results$p_val_adj == 0] <- min_nonzero_pval

# 绘制火山图
volcano_plot <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.7, size = 1.5) +
  scale_color_manual(values = c("Up in Cold" = "#D53E4F", "Down in Cold" = "#3288BD", "Not Sig" = "grey80")) +
  geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = gene), 
                  size = 4, box.padding = 0.5, point.padding = 0.3, max.overlaps = 20, color = "black") +
  theme_classic(base_size = 14) +
  labs(title = paste0("Volcano Plot: ", group_cold, " vs ", group_rt, " in Monocytes"),
       x = expression("Log"[2]*" Fold Change"),
       y = expression("-Log"[10]*" (Adjusted P-value)")) +
  theme(legend.position = "top", legend.title = element_blank(), plot.title = element_text(hjust = 0.5))

# 保存火山图
ggsave(filename = file.path("plots", "Volcano_Monocytes_Cold_vs_RT.pdf"), plot = volcano_plot, width = 8, height = 7)
ggsave(filename = file.path("plots", "Volcano_Monocytes_Cold_vs_RT.png"), plot = volcano_plot, width = 8, height = 7, dpi = 300)

# ------------------------------------------------------------------------------
# 4. 执行 GO (Biological Process) 富集分析
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在执行 GO 生物学过程 (BP) 富集分析...")

# 分别提取上调和下调的基因列表
genes_up <- deg_results %>% filter(Significance == "Up in Cold") %>% pull(gene)
genes_down <- deg_results %>% filter(Significance == "Down in Cold") %>% pull(gene)

# 定义一个运行 GO 分析的便捷函数
run_go <- function(gene_list, direction) {
  if (length(gene_list) < 5) {
    print(paste("   ⚠️", direction, "显著基因少于5个，跳过富集。"))
    return(NULL)
  }
  
  go_res <- enrichGO(gene          = gene_list,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'SYMBOL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.2,
                     readable      = FALSE)
  
  if (!is.null(go_res) && nrow(go_res@result) > 0) {
    df <- go_res@result %>% filter(p.adjust < 0.05) %>% mutate(Direction = direction)
    return(df)
  }
  return(NULL)
}

# 运行富集
go_up <- run_go(genes_up, "Up_in_Cold")
go_down <- run_go(genes_down, "Down_in_Cold")

# 合并结果并导出
all_go_results <- bind_rows(go_up, go_down)

if (nrow(all_go_results) > 0) {
  all_go_results <- all_go_results %>% dplyr::select(Direction, everything())
  output_go_file <- file.path("files", "GO_BP_Enrichment_Monocytes_Cold_vs_RT.csv")
  write.csv(all_go_results, file = output_go_file, row.names = FALSE)
  print(paste("✅ GO 富集分析完成！结果已保存至:", output_go_file))
} else {
  print("⚠️ 警告：没有找到任何显著富集的 GO 通路。")
}

print("🎉 所有分析步骤已顺利完成！")