# ==============================================================================
# Script 5: Monocytes 亚群 Marker 基因热图 (Top 20 版)
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

# 剔除 TN_30C 组（根据原脚本逻辑）
mono_cells = subset(mono_cells, Group == 'TN_30C', invert = TRUE)

# ------------------------------------------------------------------------------
# 步骤2: 筛选 Top 20 Marker 基因
# ------------------------------------------------------------------------------
print("🚀 步骤2: 筛选 Padj < 0.05 且 logFC 排名前 20 的 Marker 基因...")

top20_markers <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) # <--- 这里从 5 改成了 20

heatmap_genes <- unique(top20_markers$gene)

# ------------------------------------------------------------------------------
# 步骤3: 提取和准备表达矩阵及元数据
# ------------------------------------------------------------------------------
print("🚀 步骤3: 提取数据并处理注释标签...")

# 1. 获取 scale.data 中实际存在的基因名
all_scaled_genes <- rownames(GetAssayData(mono_cells, layer = "scale.data"))

# 2. 【关键修复】：只保留那些既在 Top 20 列表中，又在 scale.data 矩阵中的基因
genes_to_plot <- heatmap_genes[heatmap_genes %in% all_scaled_genes]

# 检查一下有多少基因被过滤掉了
missing_genes <- setdiff(heatmap_genes, genes_to_plot)
if(length(missing_genes) > 0){
  message("⚠️ 提示：有 ", length(missing_genes), " 个基因不在 scale.data 中，已被跳过。")
  # print(missing_genes) # 如果想看具体是哪些，可以取消注释
}

# 3. 提取 Scale 后的表达矩阵 (使用过滤后的基因名)
mat <- GetAssayData(mono_cells, layer = "scale.data")[genes_to_plot, ]

# 4. 提取 metadata
meta <- mono_cells@meta.data[, c("Group", "mono_cluster_id", "mono_annotation")]

# 5. 设置因子水平以保持顺序
annotation_levels <- unique(meta[order(meta$mono_cluster_id), "mono_annotation"])
meta$mono_annotation <- factor(meta$mono_annotation, levels = annotation_levels)

# 6. 对细胞进行排序
cell_order <- order(meta$mono_cluster_id, meta$Group)
mat_sorted <- mat[, cell_order]
meta_sorted <- meta[cell_order, ]


# ------------------------------------------------------------------------------
# 步骤4: 使用 ComplexHeatmap 绘制
# ------------------------------------------------------------------------------
print("🚀 步骤4: 正在生成 ComplexHeatmap (Top 20 基因)...")

# 1. 手动定义颜色映射表
my_colors <- list(
  Condition = c(
    "Cold_4C" = "#4363d8", 
    "RT_25C"  = "#3cb44b", 
    "TN_30C"  = "#e6194b"  
  ),
  Cluster = c(
    "Activated/Mature Monocytes"         = "#E41A1C",
    "Pro-inflammatory classic Monocytes" = "#377EB8",
    "FKBP5+ classic Monocytes"           = "#4DAF4A",
    "FKBP5+ non-classic Monocytes"       = "#984EA3",
    "Antigen-Presenting Monocytes"       = "#FF7F00",
    "ISG high Monocytes"                 = "#FFFF33"
  )
)

# 2. 定义顶部注释
top_anno <- HeatmapAnnotation(
  Condition = meta_sorted$Group,
  Cluster = meta_sorted$mono_annotation, 
  col = my_colors,
  show_annotation_name = TRUE,
  annotation_name_side = "right"
)

# 定义表达量颜色梯度
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
  # 【调整】：由于基因变多，将 fontsize 从 9 调小到 6-7，避免名字重叠
  row_names_gp = gpar(fontface = "italic", fontsize = 7), 
  column_split = meta_sorted$mono_cluster_id, 
  column_title = NULL,              
  border = TRUE                     
)

# ------------------------------------------------------------------------------
# 步骤5: 保存结果 (调整了图片高度以适应更多基因)
# ------------------------------------------------------------------------------
print("🚀 步骤5: 正在保存图片...")

# 【调整】：基因变多，将 height 从 9 增加到 14，确保纵向空间充裕
png_filename <- file.path("pictures", "Heatmap_Monocytes_Top20_ComplexHeatmap.png")
png(png_filename, width = 12, height = 14, units = "in", res = 300)
draw(p_heatmap, column_title = "Top 20 Marker Genes across Monocyte Subgroups", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

pdf_filename <- file.path("pictures", "Heatmap_Monocytes_Top20_ComplexHeatmap.pdf")
pdf(pdf_filename, width = 12, height = 14)
draw(p_heatmap, column_title = "Top 20 Marker Genes across Monocyte Subgroups", 
     column_title_gp = gpar(fontsize = 16, fontface = "bold"))
dev.off()

print(paste("✅ 脚本执行完毕！Top 20 热图已保存至:", png_filename))