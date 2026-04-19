# ==============================================================================
# Script 3: 提取 PBMC Monocytes、Leiden 重聚类与亚群差异分析 (整理输出版)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# ------------------------------------------------------------------------------
# 0. 环境准备与创建输出目录
# ------------------------------------------------------------------------------
library(reticulate)
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(harmony)
library(leidenbase)

# 创建输出子文件夹
if (!dir.exists("files")) dir.create("files")
if (!dir.exists("pictures")) dir.create("pictures")

# ------------------------------------------------------------------------------
# 1. 读取数据并提取 PBMC 中的 Monocytes
# ------------------------------------------------------------------------------
print("🚀 步骤1/5: 读取注释好的组织数据并提取 PBMC 中的 Monocytes...")
sc_by_tissue <- qread('sc_by_tissue.qs')

# 直接获取 PBMC 数据
pbmc_obj <- sc_by_tissue[["PBMC"]]

# 检查并提取 Monocytes
if ("Monocytes" %in% unique(pbmc_obj$SingleR.labels)) {
  mono_combined <- subset(pbmc_obj, subset = SingleR.labels == "Monocytes")
} else {
  stop("❌ 错误：在 PBMC 组织中未找到 Monocytes！")
}

# Seurat V5 子集提取后使用 JoinLayers，确保后续流程正常
mono_combined <- JoinLayers(mono_combined)

# ------------------------------------------------------------------------------
# 2. 标准降维、Harmony 批次校正与 Leiden 重聚类
# ------------------------------------------------------------------------------
print("🚀 步骤2/5: 正在进行数据归一化、降维、Harmony 批次校正与 Leiden 聚类...")

mono_combined <- NormalizeData(mono_combined, normalization.method = "LogNormalize", scale.factor = 10000)
mono_combined <- FindVariableFeatures(mono_combined, selection.method = "vst")
mono_combined <- ScaleData(mono_combined, vars.to.regress = c("nCount_RNA", "percent.mt"))

# 针对具体细胞亚群的重聚类，PC数目设为 30
mono_combined <- RunPCA(mono_combined, npcs = 30, verbose = FALSE)

# 使用 orig.ident 进行 Harmony 批次校正
mono_combined <- RunHarmony(mono_combined, group.by.vars = "orig.ident", verbose = FALSE)

mono_combined <- RunUMAP(mono_combined, reduction = "harmony", dims = 1:30)
mono_combined <- FindNeighbors(mono_combined, reduction = "harmony", dims = 1:30)

# 使用 Leiden 聚类 (algorithm = 4)
mono_combined <- FindClusters(mono_combined, resolution = 0.5, algorithm = 4)

# ------------------------------------------------------------------------------
# 3. 计算并输出差异表达 Marker (保存至 files 路径)
# ------------------------------------------------------------------------------
print("🚀 步骤3/5: 正在计算各亚群 Top 20 Markers...")
Idents(mono_combined) <- "seurat_clusters"

mono_markers <- FindAllMarkers(mono_combined, 
                               only.pos = TRUE, 
                               min.pct = 0.25, 
                               logfc.threshold = 0.25,
                               verbose = FALSE)

top20_mono_markers <- mono_markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC)

# 输出到 files 文件夹
write.csv(mono_markers, "files/Markers_All_PBMC_Monocytes_Subclusters.csv", row.names = FALSE)
write.csv(top20_mono_markers, "files/Markers_Top20_PBMC_Monocytes_Subclusters.csv", row.names = FALSE)
print("   ✅ Marker 列表已导出至 files/ 目录。")

# ------------------------------------------------------------------------------
# 4. 绘制并分别输出图片 (保存至 pictures 路径)
# ------------------------------------------------------------------------------
print("🚀 步骤4/5: 正在生成 UMAP 图谱并导出...")

# 图 1：聚类总图
p_total <- DimPlot(mono_combined, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + 
  ggtitle("PBMC Monocytes - Subclusters") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))

ggsave(filename = "pictures/PBMC_Monocytes_Total_UMAP.png", plot = p_total, width = 8, height = 7, dpi = 300)

# 图 2：Cold vs RT 对比图 (上下拼接)
p_cold <- DimPlot(subset(mono_combined, subset = Group == "Cold_4C"), reduction = "umap", group.by = "seurat_clusters", label = FALSE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  NoLegend()

p_rt <- DimPlot(subset(mono_combined, subset = Group == "RT_25C"), reduction = "umap", group.by = "seurat_clusters", label = FALSE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  NoLegend()

p_comparison <- p_cold / p_rt

ggsave(filename = "pictures/PBMC_Monocytes_Condition_Comparison_UMAP.png", plot = p_comparison, width = 7, height = 10, dpi = 300)

# ------------------------------------------------------------------------------
# 5. 保存结果对象 (保存至 files 路径)
# ------------------------------------------------------------------------------
qsave(mono_combined, "files/sc_pbmc_monocytes_reclustered.qs")
print("✅ 脚本3 运行完毕。")
print("   📁 数据结果位于: files/")
print("   🖼️ 图片结果位于: pictures/")