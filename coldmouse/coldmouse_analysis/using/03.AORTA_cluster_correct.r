# 此脚本用于纠正分群后的出图，请在02之后使用
library(Seurat)
library(qs)
library(patchwork)
library(ggplot2)

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
sc_by_tissue = qread('sc_by_tissue_louvain.qs')

# 注意：根据您上传的表格内容是Aorta（主动脉）的数据，这里将变量名更新为 aorta
# 若您原本的数据集名字不叫 Aorta，请自行修改此处的 [['Aorta']]
aorta = sc_by_tissue[['Aorta']]

# ==========================================
# 1. 剔除 cluster 7 和 16
# ==========================================
aorta <- subset(aorta, subset = seurat_clusters %in% c("7", "16"), invert = TRUE)

# ==========================================
# 2. 字典映射关系 (根据 Aorta CSV 表格)
# 注意：0, 3, 12 合并为 Smooth muscle cell; 2, 9 合并为 Fibromyocytes 等
# ==========================================
cluster2celltype <- c(
  "0" = "Smooth muscle cell",
  "1" = "proinflammation Macrophage",
  "2" = "Fibromyocytes",
  "3" = "Smooth muscle cell",
  "4" = "osteochondrogenic cell",
  "5" = "Adventitial Fibroblasts",
  "6" = "Vascular Tissue-Resident Macrophages",
  "8" = "VSMC-derived Osteoblast-like cells",
  "9" = "Fibromyocytes",
  "10" = "Lipid-Associated Macrophages",
  "11" = "Arterial Endothelial Cells",
  "12" = "Smooth muscle cell",
  "13" = "Fibroblast",
  "14" = "moDCs",
  "15" = "Cytotoxic CD8+",
  "17" = "Naive/Central Memory T",
  "18" = "Proliferating T",
  "19" = "Neutrophils",
  "20" = "Pericytes",
  "21" = "Adventitial Fibroblasts",
  "22" = "Capillary Endothelial Cells",
  "23" = "gamma delta T", # 移除了LaTeX符号以防止ggplot2报错
  "24" = "B cells"
)

# ==========================================
# 3. 为各个合并后的小类指定颜色字典
# 大类相同的用相近颜色，利用不同深浅区分小类
# ==========================================
celltype_colors_dict <- c(
  # SMC大类 (绿色)
  "Smooth muscle cell" = "#238b45",
  
  # Macrphage大类 (蓝色系，不同深浅)
  "proinflammation Macrophage" = "#c6dbef",
  "Vascular Tissue-Resident Macrophages" = "#6baed6",
  "Lipid-Associated Macrophages" = "#2171b5",
  
  # Fibromyocyte大类 (橙色)
  "Fibromyocytes" = "#fd8d3c",
  
  # VSMC-derived大类 (紫色系，不同深浅)
  "osteochondrogenic cell" = "#bcbddc",
  "VSMC-derived Osteoblast-like cells" = "#756bb1",
  
  # Fibroblasts大类 (红色系，不同深浅)
  "Adventitial Fibroblasts" = "#fcbba1",
  "Fibroblast" = "#cb181d",
  
  # Endothelial cells大类 (深绿色系，不同深浅)
  "Arterial Endothelial Cells" = "#a1d99b",
  "Capillary Endothelial Cells" = "#006d2c",
  
  # T细胞大类 (粉/紫红色系，不同深浅)
  "Cytotoxic CD8+" = "#feebe2",
  "Naive/Central Memory T" = "#fbb4b9",
  "Proliferating T" = "#f768a1",
  "gamma delta T" = "#ae017e",
  
  # 其他独立类群
  "moDCs" = "#d94801",       # 焦橙色
  "Neutrophils" = "#969696", # 灰色
  "Pericytes" = "#8c510a",   # 棕色
  "B cells" = "#01665e"      # 青色
)

# ==========================================
# 4. 合并相同 celltype 并重新顺序编号 (1, 2, 3...)
# ==========================================
# 获取基础文本注释
mapped_names <- unname(cluster2celltype[as.character(aorta$seurat_clusters)])
pure_celltypes <- ifelse(is.na(mapped_names), as.character(aorta$seurat_clusters), mapped_names)

# 提取所有不重复的细胞类型并排序，然后生成 1, 2, 3... 的新编号字典
unique_types <- sort(unique(pure_celltypes))
type2newid <- setNames(as.character(seq_along(unique_types)), unique_types)
id2type <- setNames(names(type2newid), type2newid) # 反向字典: "1" = "B cells" ...

# 将新生成的编号作为一个新列存入 aorta
aorta$new_clusters <- unname(factor(type2newid[pure_celltypes], levels = as.character(seq_along(unique_types))))

# 更新 celltype 供保存
aorta$celltype <- paste0(aorta$new_clusters, ": ", pure_celltypes)
aorta$cell_annotation = factor(pure_celltypes)
qsave(aorta, 'aorta_corrected.qs')

# 如果内存够大不需要重新 load，可以直接注释掉下一行
#aorta = qread('aorta_corrected.qs')

# ==========================================
# 5. 准备图例文本与颜色向量
# ==========================================
cluster_levels <- levels(aorta$new_clusters)

# 自动生成 "1: B cells" 这种格式的图例文字
legend_labels <- sapply(cluster_levels, function(x) {
  paste0(x, ": ", id2type[x])
})

# 按当前 cluster 1, 2, 3... 的顺序，提取匹配的颜色，以确保颜色映射准确
final_colors <- unname(celltype_colors_dict[id2type[cluster_levels]])

# ==========================================
# 6. 绘图 (基于 scale_color_manual 使用自定义色系)
# ==========================================
p_total_anno <- DimPlot(aorta, reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_manual(values = final_colors, labels = legend_labels) + 
  ggtitle('Corrected Aorta') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_anno <- DimPlot(subset(aorta, subset = Group == "Cold_4C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_manual(values = final_colors) +
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_anno <- DimPlot(subset(aorta, subset = Group == "RT_25C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_manual(values = final_colors) +
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_anno <- DimPlot(subset(aorta, subset = Group == "TN_30C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_manual(values = final_colors) +
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 合并图层
p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

# 保存图片 (文件名更新为 Aorta 避免和旧图片混淆)
if(!dir.exists("pictures")) dir.create("pictures")
file_name_anno_png <- file.path("pictures", "UMAP_Grid_corrected_Aorta.png")
ggsave(filename = file_name_anno_png, plot = p_final_anno, width = 15, height = 11, dpi = 300)

file_name_anno_pdf <- file.path("pictures", "UMAP_Grid_corrected_Aorta_update1.pdf")
ggsave(filename = file_name_anno_pdf, plot = p_final_anno, width = 15, height = 11, dpi = 300)
# ==========================================
# 6. 生成 Cluster x Marker 气泡图 (Top 5 Markers)
# ==========================================
# 加载处理数据框所需的包
library(dplyr)

# 设置当前分类标签为你刚刚生成的合并注释列 'celltype'
Idents(aorta) <- "cell_annotation"

# 计算所有 cluster 的 marker 基因 
# (注意：如果数据较大会运行几分钟，only.pos = TRUE 只找高表达基因，加速计算)
message("Finding all markers for each cluster...")
aorta.markers <- FindAllMarkers(aorta, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# 提取每个 cluster 按照 log2FC 排序的前五的 marker 基因
top5_markers <- aorta.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

# 提取去重后的基因名列表，保证 X 轴的基因顺序与 cluster 排序大致对应
top5_genes <- unique(top5_markers$gene)

# ==========================================
# 绘制 DotPlot (气泡图)
# ==========================================
p_dotplot <- DotPlot(aorta, features = top5_genes, group.by = "cell_annotation") +
  # 调整颜色渐变，模仿附图中的由白到深红
  scale_color_gradient(low = "white", high = "darkred") + 
  # 修改图例标题以贴合附图
  labs(
    size = "Fraction of cells\nin group (%)", 
    color = "Mean expression\nin group"
  ) +
  theme_bw() + # 基础黑白主题去灰底
  ggtitle('AORTA Dimplot')+
  theme(
    # X轴基因名字体倾斜 45 度，防止基因名重叠
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "italic", color = "black"),
    axis.text.y = element_text(color = "black"),
    # 移除横纵轴默认的 "Features" 和 "Identity" 标题
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    # 调整图例背景和边框
    panel.grid.major = element_line(color = "grey90", size = 0.5),
    panel.grid.minor = element_blank()
  )

# ==========================================
# 保存图片
# ==========================================
file_name_dotplot_png <- file.path("pictures", "DotPlot_AORTA_Top5_Markers.png")
ggsave(filename = file_name_dotplot_png, plot = p_dotplot, width = 20, height = 6, dpi = 300)

# file_name_dotplot_pdf <- file.path("pictures", "DotPlot_Top5_Markers.pdf")
# ggsave(filename = file_name_dotplot_pdf, plot = p_dotplot, width = 16, height = 6, dpi = 300)

message("DotPlot saved to pictures directory.")