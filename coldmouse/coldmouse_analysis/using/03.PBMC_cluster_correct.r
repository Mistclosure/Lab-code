#此脚本用于纠正分群后的出图，请在02之后使用
library(Seurat)
library(qs)
library(patchwork)
library(ggplot2)

setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
sc_by_tissue = qread('sc_by_tissue_louvain.qs')
pbmc = sc_by_tissue[['PBMC']]

# 1. 剔除 cluster 26 和 27
pbmc <- subset(pbmc, subset = seurat_clusters %in% c("26", "27"), invert = TRUE)

# 2. 字典映射关系
cluster2celltype <- c(
  "0" = "Mature B cells", "1" = "NK cells", "2" = "Monocytes", 
  "3" = "Monocytes", "4" = "Naive CD4", "5" = "Effector CD8 cells", 
  "6" = "Monocytes", "7" = "Platelets", "8" = "Tregs", 
  "9" = "Naive CD8", "10" = "Mature B cells", "11" = "Effector CD8 cells", 
  "12" = "Monocytes", "13" = "Basophils", "14" = "Mature B cells", 
  "15" = "Exhausted CD8", "16" = "Naive Bcells", "17" = "Resting T cells", 
  "18" = "Tregs", "19" = "Erythrocytes", "20" = "Infiltrating Neutrophils", 
  "21" = "Dendritic cell", "22" = "Monocytes", "23" = "gamma delta", 
  "24" = "Activated Neutrophils", "25" = "gamma delta T17 cells", 
  "28" = "Epithelial cells"
)

# ==========================================
# 3. 合并相同 celltype 并重新顺序编号
# ==========================================
# 先获取基础文本注释（没匹配到的保留原 cluster 编号）
mapped_names <- unname(cluster2celltype[as.character(pbmc$seurat_clusters)])
pure_celltypes <- ifelse(is.na(mapped_names), as.character(pbmc$seurat_clusters), mapped_names)

# 提取所有不重复的细胞类型并排序，然后生成 1, 2, 3... 的新编号字典
unique_types <- sort(unique(pure_celltypes))
type2newid <- setNames(as.character(seq_along(unique_types)), unique_types)
id2type <- setNames(names(type2newid), type2newid) # 反向字典，供稍后生成图例使用

# 将新生成的 1, 2, 3... 编号作为一个新列存入 pbmc
# 这里使用 factor 确保画图时顺序严格按照 1, 2, 3 排列（防止 1 后接 10）
pbmc$new_clusters <- unname(factor(type2newid[pure_celltypes], levels = as.character(1:length(unique_types))))

# 同样更新 celltype 供保存
pbmc$celltype <- paste0(pbmc$new_clusters, ": ", pure_celltypes)
qsave(pbmc, 'pbmc_corrected.qs')
pbmc = qread('pbmc_corrected.qs')
# ==========================================
# 4. 进一步过滤和合并群 (基于旧的 new_clusters)
# ==========================================
# A. 剔除指定群
pbmc <- subset(pbmc, subset = new_clusters %in% c("1", "5", "6", "9", "13", "18"), invert = TRUE)

# B. 合并 4 和 15 为 CD8 (此时仍使用旧编号)
pbmc$new_clusters <- as.character(pbmc$new_clusters)
pbmc$new_clusters[pbmc$new_clusters == "15"] <- "4"
id2type["4"] <- "CD8" # 确保 4 号的定义更新为 CD8

# ==========================================
# 5. 核心：重新计算连续序号
# ==========================================
# 获取当前 pbmc 中实际存在的旧编号，并按数字大小排序
old_ids_present <- sort(as.numeric(unique(pbmc$new_clusters)))
old_ids_present <- as.character(old_ids_present)

# 获取这些旧编号对应的细胞类型名称
current_types <- id2type[old_ids_present]

# 创建重映射字典：旧编号 -> 新连续编号 (1, 2, 3...)
new_id_map <- setNames(as.character(seq_along(old_ids_present)), old_ids_present)

# 应用重映射到 Seurat 对象
pbmc$new_clusters <- unname(factor(new_id_map[pbmc$new_clusters], 
                             levels = as.character(seq_along(old_ids_present))))

# 更新 id2type 映射表，供后续 legend_labels 使用
# 现在 id2type 变成了：新连续编号 -> 细胞类型
id2type <- setNames(current_types, as.character(seq_along(old_ids_present)))

# 同步更新 celltype 字符串列
pbmc$celltype <- paste0(pbmc$new_clusters, ": ", id2type[as.character(pbmc$new_clusters)])

# ==========================================
# 准备替换图例文本的向量
# ==========================================
cluster_levels <- levels(pbmc$new_clusters)

# 自动生成 "1: Monocytes" 这种格式的图例文字
legend_labels <- sapply(cluster_levels, function(x) {
  paste0(x, ": ", id2type[x])
})

# ==========================================
# 绘图: 基于合并后的 new_clusters 绘图
# ==========================================
# 将 group.by 改为 "new_clusters"，图上打出的 Label 就会变成全新的 1,2,3...
p_total_anno <- DimPlot(pbmc, reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  scale_color_discrete(labels = legend_labels) + # 替换图例的文字显示
  ggtitle('Corrected PBMC') + 
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())

p_cold_anno <- DimPlot(subset(pbmc, subset = Group == "Cold_4C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  ggtitle("Cold_4C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_rt_anno <- DimPlot(subset(pbmc, subset = Group == "RT_25C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  ggtitle("RT_25C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

p_tn_anno <- DimPlot(subset(pbmc, subset = Group == "TN_30C"), reduction = "umap", group.by = "new_clusters", label = TRUE, repel = TRUE) + 
  ggtitle("TN_30C") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
  NoLegend()

# 合并图层
p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
  plot_layout(guides = "collect", axes = "collect") & 
  theme(legend.text = element_text(size = 9))

if(!dir.exists("pictures")) dir.create("pictures")
file_name_anno_png <- file.path("pictures", "UMAP_Grid_corrected_PBMC.png")
ggsave(filename = file_name_anno_png, plot = p_final_anno, width = 15, height = 11, dpi = 300)
file_name_anno_pdf <- file.path("pictures", "UMAP_Grid_corrected_PBMC_update1.pdf")
ggsave(filename = file_name_anno_pdf, plot = p_final_anno, width = 15, height = 11, dpi = 300)
