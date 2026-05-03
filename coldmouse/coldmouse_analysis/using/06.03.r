# ==============================================================================
# Script 06.03: Intra-sample Trajectory (PBMC Monocytes Only)
# Purpose: Research Monocyte state transitions within PBMC (e.g. RT -> Cold)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
library(qs)
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)
library(patchwork) 
library(scales)    

# ------------------------------------------------------------------------------
# 【参数配置区】
# ------------------------------------------------------------------------------
target_gene <- "Fkbp5"       # 🌟 需要额外标注的基因名 (请确认是否为 Fkbp5 等)
highlight_color <- "red"    # 🌟 该基因表达量 > 0 的细胞的高亮颜色

# 1. 读取数据并提取 PBMC Monocytes
print("Loading data and subsetting PBMC Monocytes...")
pbmc <- qread('pbmc_corrected.qs')
pbmc_mono <- subset(pbmc, new_clusters == "12") 

# 2. 构建 Monocle 3 对象
print("Converting to Monocle 3...")
expression_matrix <- LayerData(pbmc_mono, layer = "counts")
cell_metadata <- pbmc_mono@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                            row.names = rownames(expression_matrix))

cds_mono <- new_cell_data_set(expression_matrix,
                              cell_metadata = cell_metadata,
                              gene_metadata = gene_metadata)
cds_mono <- estimate_size_factors(cds_mono)

# 3. 降维与聚类
print("Preprocessing and dimension reduction...")
cds_mono <- preprocess_cds(cds_mono, num_dim = 30)

print("Aligning Group differences...")
cds_mono <- align_cds(cds_mono, alignment_group = "Group")

print("Reducing dimensions...")
cds_mono <- reduce_dimension(cds_mono, 
                             reduction_method = "UMAP", 
                             preprocess_method = "Aligned", 
                             umap.n_neighbors = 50,         
                             umap.min_dist = 0.5)           

cds_mono <- cluster_cells(cds_mono, reduction_method = "UMAP")
cds_mono <- learn_graph(cds_mono, use_partition = FALSE)

# 4. 学习轨迹图
print("Learning Trajectory Graph...")
cds_mono <- learn_graph(cds_mono, use_partition = FALSE, 
                        learn_graph_control = list(minimal_branch_len = 5))

# 5. 设定拟时序起点并排序
print("Ordering cells based on Early Marker (Ly6c2)...")
get_earliest_principal_node <- function(cds, marker_gene = "Ly6c2"){
  gene_expr <- normalized_counts(cds)[marker_gene, ]
  threshold <- quantile(gene_expr, 0.95)
  early_cells <- names(gene_expr[gene_expr >= threshold])
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  early_nodes <- as.numeric(names(which.max(table(closest_vertex[early_cells, ]))))
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[early_nodes]
  return(root_pr_nodes)
}

root_node <- get_earliest_principal_node(cds_mono, marker_gene = "Ly6c2")
cds_mono <- order_cells(cds_mono, root_pr_nodes = root_node)

# ------------------------------------------------------------------------------
# 6. 可视化 (🌟 基因亚群高亮逻辑)
# ------------------------------------------------------------------------------
print("Generating Plots...")

# 提取目标基因表达量 > 0 的细胞条形码
if (target_gene %in% rownames(cds_mono)) {
  gene_expr <- normalized_counts(cds_mono)[target_gene, ]
  target_pos_cells <- names(which(gene_expr > 0))
  print(paste("Found", length(target_pos_cells), "cells expressing", target_gene))
} else {
  warning(paste("⚠️ 基因", target_gene, "在数据中未找到，请检查拼写。高亮将失效。"))
  target_pos_cells <- character(0)
}

# (A) 构建高亮拼图
all_groups <- as.character(unique(colData(cds_mono)$Group))
base_colors <- hue_pal()(length(all_groups))
color_mapping <- setNames(base_colors, all_groups)
color_mapping["Other"] <- "#d3d3d3" 

# 🌟 为目标基因阳性的细胞新增一种颜色映射
target_label <- paste0("Cold_4C (", target_gene, "+)")
color_mapping[target_label] <- highlight_color 

plot_list <- list()

for (grp in all_groups) {
  highlight_col_name <- paste0("HL_", grp)
  
  # 默认分配：是当前组就填组名，否则填 Other
  cell_labels <- ifelse(colData(cds_mono)$Group == grp, grp, "Other")
  
  # 🌟 针对 Cold_4C 组的特殊判定逻辑
  if (grp == "Cold_4C") {
    # 找到既属于 Cold_4C，又表达该基因的细胞
    is_target_pos <- (colData(cds_mono)$Group == "Cold_4C") & (colnames(cds_mono) %in% target_pos_cells)
    cell_labels[is_target_pos] <- target_label
    
    # 设置因子水平：Other垫底 -> 普通Cold居中 -> 阳性Cold放顶层
    factor_levels <- c("Other", grp, target_label)
  } else {
    factor_levels <- c("Other", grp)
  }
  
  colData(cds_mono)[[highlight_col_name]] <- factor(cell_labels, levels = factor_levels)
  
  # 重新排序确保目标细胞最后绘制（浮在最上面）
  cds_ordered <- cds_mono[, order(colData(cds_mono)[[highlight_col_name]])]
  
  p <- plot_cells(cds_ordered, color_cells_by = highlight_col_name, 
                  label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
                  cell_size = 0.8) + 
    scale_color_manual(values = color_mapping) +
    ggtitle(grp) +
    # 保留图例以便认出红色点代表什么，但不显示图例标题
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
          legend.title = element_blank())
  
  plot_list[[grp]] <- p
}

# 拼合图片，并将图例统一收集到下方
plot_group_highlight <- wrap_plots(plot_list, ncol = length(all_groups)) + 
  plot_layout(guides = "collect") +
  plot_annotation(title = "PBMC Monocytes: Group Trajectory Highlights") &
  theme(legend.position = "bottom")

# (B) 其他图保持原样
plot_time <- plot_cells(cds_mono, color_cells_by = "pseudotime", 
                        label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("PBMC Monocytes: Pseudotime")

plot_split <- plot_cells(cds_mono, color_cells_by = "pseudotime", 
                         label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  facet_wrap(~Group) + 
  ggtitle("PBMC Monocytes: Pseudotime Split by Group")

# 7. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

#qsave(cds_mono, "files/script_06_03_monocle3_pbmc_monocytes.qs")

ggsave("pictures/monocle3_06_03_group_highlight.pdf", plot_group_highlight, width = 4 * length(all_groups), height = 5.5)
ggsave("pictures/monocle3_06_03_pseudotime.pdf", plot_time, width = 6, height = 5)
ggsave("pictures/monocle3_06_03_split.pdf", plot_split, width = 10, height = 5)

print("Script 06.03 COMPLETED successfully!")