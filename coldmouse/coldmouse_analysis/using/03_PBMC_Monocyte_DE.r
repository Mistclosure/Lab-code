# ==============================================================================
# Script 3: PBMC Monocyte 亚群差异表达分析 (DEA)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)

# 1. 加载数据 (默认加载 Louvain 版本，可根据需要修改)
cluster_method <- "louvain" 
input_file <- paste0('sc_by_tissue_', cluster_method, '.qs')

if (!file.exists(input_file)) {
  stop(paste("未找到输入文件:", input_file, "请先确保 Script 2 运行成功。"))
}

sc_by_tissue <- qread(input_file)

# 2. 提取 PBMC 对象
pbmc <- sc_by_tissue[["PBMC"]]

# 3. 提取 Monocyte 亚群
# 注意：SingleR 的标签通常是 "Monocytes"，请确保名称与输出的 Marker 列表一致
monocyte_label <- "Monocytes" # 如果你的标签是 "Monocyte" 请在此修改

if (!(monocyte_label %in% unique(pbmc$SingleR.labels))) {
  print("当前的标签列表包含：")
  print(unique(pbmc$SingleR.labels))
  stop(paste("在 PBMC 中未找到标签:", monocyte_label))
}

pbmc_mono <- subset(pbmc, subset = SingleR.labels == monocyte_label)

# 切换标识符为分组信息 (Group)
Idents(pbmc_mono) <- "Group"

print(paste("🚀 开始对 PBMC", monocyte_label, "进行差异分析..."))

# ------------------------------------------------------------------------------
# 4. 差异分析：Cold_4C vs RT_25C
# ------------------------------------------------------------------------------
print(">>> 正在计算: Cold_4C vs RT_25C ...")
deg_cold_vs_rt <- FindMarkers(
  pbmc_mono, 
  ident.1 = "Cold_4C", 
  ident.2 = "RT_25C",
  test.use = "wilcox", # 默认使用 Wilcoxon 秩和检验
  min.pct = 0.1,       # 至少在 10% 的细胞中表达
  logfc.threshold = 0.25
)

deg_cold_vs_rt$gene <- rownames(deg_cold_vs_rt)
write.csv(deg_cold_vs_rt, 
          file.path("files", paste0("DEG_PBMC_Mono_Cold_vs_RT_", cluster_method, ".csv")), 
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 5. 差异分析：Cold_4C vs TN_30C
# ------------------------------------------------------------------------------
print(">>> 正在计算: Cold_4C vs TN_30C ...")
deg_cold_vs_tn <- FindMarkers(
  pbmc_mono, 
  ident.1 = "Cold_4C", 
  ident.2 = "TN_30C",
  test.use = "wilcox",
  min.pct = 0.1,
  logfc.threshold = 0.25
)

deg_cold_vs_tn$gene <- rownames(deg_cold_vs_tn)
write.csv(deg_cold_vs_tn, 
          file.path("files", paste0("DEG_PBMC_Mono_Cold_vs_TN_", cluster_method, ".csv")), 
          row.names = FALSE)

print("✅ 差异分析完成，结果已保存至 files 文件夹。")

# ------------------------------------------------------------------------------
# 6. (可选) 简单可视化：查看核心差异基因
# ------------------------------------------------------------------------------
# 打印 Cold vs RT 前 10 个上调基因
print("Top 10 Up-regulated genes in Cold vs RT:")
print(head(deg_cold_vs_rt %>% filter(avg_log2FC > 0) %>% arrange(p_val_adj), 10))