# ==============================================================================
# Script 5: PBMC Cold_4C 单核细胞组内异质性分析 (Fkbp5+ vs Fkbp5-)
# 独立运行，直接依赖数据: sc_by_tissue.qs (脚本2产物)
# 包含：差异基因分析 (DEA) -> 火山图绘制 -> GO/KEGG 通路富集分析
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠基因数据库
library(enrichplot)

# ------------------------------------------------------------------------------
# 1. 读取数据并严格提取目标细胞亚群 (基于全局 SingleR 注释)
# ------------------------------------------------------------------------------
print("🚀 正在加载 sc_by_tissue.qs ...")
sc_by_tissue <- qread("sc_by_tissue.qs")

# 检查 PBMC 是否存在
if (!"PBMC" %in% names(sc_by_tissue)) stop("❌ 未找到 PBMC 数据！")

pbmc_obj <- sc_by_tissue[["PBMC"]]
pbmc_obj <- JoinLayers(pbmc_obj) # 确保 Seurat V5 矩阵层已合并

# 提取条件：Group == Cold_4C 且 SingleR.labels == "Monocytes"
print("🚀 正在提取 PBMC 中 Cold_4C 组的单核细胞...")
target_cells <- subset(pbmc_obj, subset = Group == "Cold_4C" & SingleR.labels == "Monocytes")

if (ncol(target_cells) < 10) stop("❌ 符合条件的细胞数太少，无法进行差异分析！")

# ------------------------------------------------------------------------------
# 2. 定义 Fkbp5 阴阳性分组 (+ / -)
# ------------------------------------------------------------------------------
# 确定靶基因名字
target_gene <- grep("Fkbp5", rownames(target_cells), ignore.case = TRUE, value = TRUE)[1]

# 提取 counts 表达量，大于 0 为阳性
fkbp5_counts <- GetAssayData(target_cells, assay = "RNA", layer = "counts")[target_gene, ]

# 【关键修改】直接使用 Fkbp5+ 和 Fkbp5- 作为内部 Metadata 标签
target_cells$Fkbp5_Status <- ifelse(fkbp5_counts > 0, "Fkbp5+", "Fkbp5-")

print("✅ 分组完成，细胞数量统计：")
print(table(target_cells$Fkbp5_Status))

# 设置为当前计算标识 (Idents)
Idents(target_cells) <- "Fkbp5_Status"

# ------------------------------------------------------------------------------
# 3. 差异表达分析 (Fkbp5+ vs Fkbp5-)
# ------------------------------------------------------------------------------
print("🚀 正在计算差异表达基因 (Fkbp5+ vs Fkbp5-)...")
# 默认正值表示在 Fkbp5+ 中上调，负值表示在 Fkbp5- 中上调
deg_results <- FindMarkers(
  target_cells, 
  ident.1 = "Fkbp5+", 
  ident.2 = "Fkbp5-", 
  logfc.threshold = 0.25, # log2FC 阈值
  min.pct = 0.05           # 至少在 10% 的细胞中表达
)

# 将基因名变为一列并保存
deg_results <- deg_results %>% 
  rownames_to_column(var = "gene") %>%
  arrange(desc(avg_log2FC))

write.csv(deg_results, "PBMC_Cold4C_Monocytes_Fkbp5Plus_vs_Minus_DEGs.csv", row.names = FALSE)
print("✅ 差异基因计算完成，已导出 CSV！")

# --- 绘制定制版火山图 (基于 ggplot2) ---
deg_results$Significance <- "Not Sig"
deg_results$Significance[deg_results$avg_log2FC > 0.25 & deg_results$p_val_adj < 0.05] <- "Up in Fkbp5+"
deg_results$Significance[deg_results$avg_log2FC < -0.25 & deg_results$p_val_adj < 0.05] <- "Up in Fkbp5-"

p_volcano <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
  geom_point(alpha = 0.8, size = 1.5) +
  scale_color_manual(values = c("Not Sig" = "grey", "Up in Fkbp5+" = "#E64B35", "Up in Fkbp5-" = "#4DBBD5")) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  theme_classic() +
  ggtitle("Fkbp5+ vs Fkbp5- Monocytes in Cold_4C PBMC") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("PBMC_Cold4C_Monocytes_Fkbp5_Volcano.png", plot = p_volcano, width = 7, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. 功能通路富集分析 (GO - 生物学过程 BP)
# ------------------------------------------------------------------------------
print("🚀 正在进行 GO (BP) 富集分析...")

# 提取在 Fkbp5+ 细胞中显著上调的基因群
genes_up <- deg_results %>% 
  filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>% 
  pull(gene)

if (length(genes_up) > 10) {
  # 运行 GO 富集
  go_bp_up <- enrichGO(
    gene          = genes_up,
    OrgDb         = org.Mm.eg.db,
    keyType       = 'SYMBOL',
    ont           = "BP", # BP = Biological Process
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2
  )
  
  if (!is.null(go_bp_up) && nrow(go_bp_up@result) > 0) {
    p_go_up <- dotplot(go_bp_up, showCategory = 15, title = "GO BP Enriched in Fkbp5+ Monocytes")
    ggsave("PBMC_Cold4C_Monocytes_Fkbp5Plus_GO_BP_Up.png", plot = p_go_up, width = 8, height = 7, dpi = 300)
    write.csv(as.data.frame(go_bp_up), "PBMC_Cold4C_Monocytes_Fkbp5Plus_GO_BP_Up.csv", row.names = FALSE)
    print("   ✅ Fkbp5+ 上调基因 GO 富集绘图完成！")
  } else {
    print("   ⚠️ Fkbp5+ 上调基因未富集到显著的 GO 通路。")
  }
} else {
  print("   ⚠️ Fkbp5+ 上调基因太少 (<10)，跳过富集分析。")
}

print("🎉 脚本 5 运行完毕！差异表、火山图和通路图已全部保存。")