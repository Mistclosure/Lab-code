# ==============================================================================
# Script 5: Monocytes 亚群 Marker 基因热图 (ComplexHeatmap 出版级复刻)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(qs)
library(ComplexHeatmap) # 核心画图包
library(circlize)       # 用于设置热图颜色映射

print("🚀 步骤1: 加载单核细胞数据与 Marker 列表...")
mono_cells <- qread("pbmc_monocytes_sub-clustered.qs")
markers_all <- read.csv(file.path("files", "Markers_All_Leiden_Monocytes_Subgroups.csv"))
mono_cells = subset(mono_cells,Group == 'TN_30C',invert = TRUE)
# ------------------------------------------------------------------------------
# 步骤2: 筛选 Top 5 Marker 基因
# ------------------------------------------------------------------------------
print("🚀 步骤2: 筛选 Padj < 0.05 且 logFC 排名前五的 Marker 基因...")

top5_markers <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 5)

heatmap_genes <- unique(top5_markers$gene)

# ------------------------------------------------------------------------------
# 步骤3: 提取和准备表达矩阵及元数据
# ------------------------------------------------------------------------------
print("🚀 步骤3: 提取数据并处理注释标签...")

# 1. 提取 Scale 后的表达矩阵
mat <- GetAssayData(mono_cells, layer = "scale.data")[heatmap_genes, ]

# 2. 提取 metadata，这次多拿一列 mono_annotation
meta <- mono_cells@meta.data[, c("Group", "mono_cluster_id", "mono_annotation")]

# 3. 【关键改动】：将 mono_annotation 设为因子，防止图例顺序错乱
# 我们利用 mono_cluster_id 的顺序来定义 mono_annotation 的级别
annotation_levels <- unique(meta[order(meta$mono_cluster_id), "mono_annotation"])
meta$mono_annotation <- factor(meta$mono_annotation, levels = annotation_levels)

# 4. 对细胞进行排序 (依然按 ID 排序，保证热图色块连续)
cell_order <- order(meta$mono_cluster_id, meta$Group)
mat_sorted <- mat[, cell_order]
meta_sorted <- meta[cell_order, ]


# ------------------------------------------------------------------------------
# 步骤4: 使用 ComplexHeatmap 绘制
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在生成带有自定义颜色的 ComplexHeatmap...")

# 1. 手动定义颜色映射表 (使用十六进制颜色码)
# 你可以根据自己的喜好修改这里的颜色
my_colors <- list(
  # 定义实验组的颜色
  Condition = c(
    "Cold_4C" = "#4363d8", # 蓝色代表冷
    "RT_25C"  = "#3cb44b", # 绿色代表常温
    "TN_30C"  = "#e6194b"  # 红色代表热
  ),
  # 定义各个单核细胞亚群的颜色 (注意：键名必须和 mono_annotation 里的文字一模一样)
  Cluster = c(
    "Activated/Mature Monocytes"         = "#E41A1C",
    "Pro-inflammatory classic Monocytes" = "#377EB8",
    "FKBP5+ classic Monocytes"           = "#4DAF4A",
    "FKBP5+ non-classic Monocytes"       = "#984EA3",
    "Antigen-Presenting Monocytes"       = "#FF7F00",
    "ISG high Monocytes"                 = "#FFFF33"
  )
)

# 2. 将定义好的颜色传入 top_anno
top_anno <- HeatmapAnnotation(
  Condition = meta_sorted$Group,
  Cluster = meta_sorted$mono_annotation, 
  col = my_colors,                  # <--- 【新增】：在这里传入你的颜色列表
  show_annotation_name = TRUE,
  annotation_name_side = "right"
)

col_fun <- colorRamp2(c(-2, 0, 2), c("#1E90FF", "white", "#B22222"))

p_heatmap <- Heatmap(
  mat_sorted,
  name = "Expression",              
  col = col_fun,                    
  top_annotation = top_anno,        
  cluster_rows = FALSE,             
  cluster_columns = FALSE,          
  show_column_names = FALSE,        
  row_names_side = "right",         
  row_names_gp = gpar(fontface = "italic", fontsize = 9), 
  # 保持用 ID 进行物理切分，这样比较整齐，且不会因为文字长短影响切片
  column_split = meta_sorted$mono_cluster_id, 
  column_title = NULL,              
  border = TRUE                     
)

# ------------------------------------------------------------------------------
# 步骤5: 保存结果
# ------------------------------------------------------------------------------
print("🚀 步骤5: 正在保存图片...")

# ComplexHeatmap 需要用 pdf() 或 png() 等设备函数配合 draw() 保存
png_filename <- file.path("pictures", "Heatmap_Monocytes_ComplexHeatmap.png")
png(png_filename, width = 12, height = 9, units = "in", res = 300)
draw(p_heatmap, column_title = "Top 5 Marker Genes across Monocyte Subgroups", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

pdf_filename <- file.path("pictures", "Heatmap_Monocytes_ComplexHeatmap.pdf")
pdf(pdf_filename, width = 12, height = 9)
draw(p_heatmap, column_title = "Top 5 Marker Genes across Monocyte Subgroups", column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

print(paste("✅ 脚本 05 执行完毕！完美复刻热图已保存至:", png_filename))