# ==============================================================================
# Script 06.02: Intra-sample Trajectory (Aorta Macrophages Only)
# Purpose: Research Macrophage state transitions within Aorta (e.g. RT -> Cold)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
library(qs)
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)

# 1. 读取数据并提取 Aorta Macrophages
print("Loading data and subsetting Aorta Macrophages...")
sc_by_tissue <- qread('sc_by_tissue_louvain.qs')
aorta <- sc_by_tissue[['Aorta']]
#aorta = subset(aorta, subset = aorta$Group == 'Cold_4C' )
aorta_mac <- subset(aorta, SingleR.labels %in% c("Macrophages"))

# 2. 构建 Monocle 3 对象
print("Converting to Monocle 3...")
expression_matrix <- LayerData(aorta_mac, layer = "counts")
cell_metadata <- aorta_mac@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                            row.names = rownames(expression_matrix))

cds_mac <- new_cell_data_set(expression_matrix,
                             cell_metadata = cell_metadata,
                             gene_metadata = gene_metadata)
cds_mac <- estimate_size_factors(cds_mac)

# 3. 降维与聚类
print("Preprocessing and dimension reduction...")
cds_mac <- preprocess_cds(cds_mac, num_dim = 30)
cds_mac <- reduce_dimension(cds_mac, reduction_method = "UMAP")
cds_mac <- cluster_cells(cds_mac, reduction_method = "UMAP", resolution = 1e-3)

# 4. 学习轨迹图
print("Learning Trajectory Graph...")
cds_mac <- learn_graph(cds_mac, use_partition = FALSE, 
                       learn_graph_control = list(minimal_branch_len = 5))

# 5. 设定拟时序起点并排序 (基于早期基因 Ly6c2)
print("Ordering cells based on Early Marker (Ly6c2)...")

get_earliest_principal_node <- function(cds, marker_gene = "Ly6c2"){
  # 提取标准化后的表达量
  gene_expr <- normalized_counts(cds)[marker_gene, ]
  # 找到表达量最高的前 5% 的细胞
  threshold <- quantile(gene_expr, 0.95)
  early_cells <- names(gene_expr[gene_expr >= threshold])
  
  # 寻找这些细胞在 UMAP 空间中最接近的轨迹节点 (Principal node)
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  
  # 统计这些高表达细胞最密集的节点序号
  early_nodes <- as.numeric(names(which.max(table(closest_vertex[early_cells, ]))))
  
  # 获取该节点的名称
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[early_nodes]
  return(root_pr_nodes)
}

# 计算起点节点并排序
root_node <- get_earliest_principal_node(cds_mac, marker_gene = "Ly6c2")
cds_mac <- order_cells(cds_mac, root_pr_nodes = root_node)

# 6. 可视化 (不显示序号，保持清爽)
print("Generating Plots...")

common_args <- list(
  cds = cds_mac,
  label_cell_groups = FALSE,       # 关闭分组标签
  label_leaves = FALSE,            # 关闭叶节点序号
  label_branch_points = FALSE,     # 关闭分支点序号
  label_roots = TRUE            # 关闭起点序号
  #label_principal_points = FALSE,  # 关闭所有主图节点序号
  #graph_label_size = 0             # 强制标签大小为0
)

plot_group <- do.call(plot_cells, c(common_args, list(color_cells_by = "Group"))) + 
  ggtitle("Aorta Macrophages: Trajectory by Group")

plot_time <- do.call(plot_cells, c(common_args, list(color_cells_by = "pseudotime"))) + 
  ggtitle("Aorta Macrophages: Pseudotime (Root: Ly6c2 high)")

plot_split <- do.call(plot_cells, c(common_args, list(color_cells_by = "pseudotime"))) + 
  facet_wrap(~Group) + 
  ggtitle("Aorta Macrophages: Pseudotime Split by Group")

# 7. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

ggsave("pictures/monocle3_aorta_macrophages_group.png", plot_group, width = 6, height = 5)
ggsave("pictures/monocle3_aorta_macrophages_pseudotime.png", plot_time, width = 6, height = 5)
#ggsave("pictures/monocle3_aorta_macrophages_split.png", plot_split, width = 10, height = 5)

print("Script 06.02 COMPLETED successfully!")