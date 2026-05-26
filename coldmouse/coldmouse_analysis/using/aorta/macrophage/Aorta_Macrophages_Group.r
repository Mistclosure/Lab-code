# ==============================================================================
# Script: Group-wise Trajectory Analysis (Aorta Macrophages Only)
# Purpose: 分组独立推断主动脉巨噬细胞的拟时序轨迹
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
library(qs)
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)

# 1. 读取数据并提取 Aorta Macrophages
print("Loading data and subsetting Aorta Macrophages...")
aorta <- qread('aorta_corrected.qs')
aorta_mac <- subset(aorta, cell_annotation %in% c("Lipid-Associated Macrophages",'proinflammation Macrophage','Vascular Tissue-Resident Macrophages'))

# 创建输出目录
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

# 2. 定义寻找起点的辅助函数 (基于早期基因 Ly6c2)
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

# 3. 设定绘图的通用参数（去掉了绑定的 cds 对象，改在循环中传入）
common_args <- list(
  label_cell_groups = FALSE,       # 关闭分组标签
  label_leaves = FALSE,            # 关闭叶节点序号
  label_branch_points = FALSE,     # 关闭分支点序号
  label_roots = TRUE               # 显示起点序号
)

# 4. 按 Group 分别进行拟时序分析
print("Running individual trajectory analysis for each Group...")
groups <- unique(aorta_mac$Group)

for (grp in groups) {
  print(paste0(">>> Processing Group: ", grp, " <<<"))
  
  # 4.1 提取当前分组的细胞
  aorta_mac_sub <- subset(aorta_mac, subset = Group == grp)
  
  # 4.2 构建单组的 Monocle 3 对象
  expr_matrix_sub <- LayerData(aorta_mac_sub, layer = "counts")
  cell_meta_sub <- aorta_mac_sub@meta.data
  gene_meta_sub <- data.frame(gene_short_name = rownames(expr_matrix_sub), 
                              row.names = rownames(expr_matrix_sub))
  
  cds_sub <- new_cell_data_set(expr_matrix_sub,
                               cell_metadata = cell_meta_sub,
                               gene_metadata = gene_meta_sub)
  cds_sub <- estimate_size_factors(cds_sub)
  
  # 4.3 降维与聚类
  # 注意：如果某个group细胞数极少，num_dim 可能需要调小，这里默认沿用 30
  cds_sub <- preprocess_cds(cds_sub, num_dim = 30)
  cds_sub <- reduce_dimension(cds_sub, reduction_method = "UMAP")
  cds_sub <- cluster_cells(cds_sub, reduction_method = "UMAP", resolution = 1e-3)
  
  # 4.4 学习轨迹图
  cds_sub <- learn_graph(cds_sub, use_partition = FALSE, 
                         learn_graph_control = list(minimal_branch_len = 5))
  
  # 4.5 设定拟时序起点并排序
  # 使用 tryCatch 处理可能出现的异常 (如某一组中几乎没有 Ly6c2 表达导致报错)
  root_node_sub <- tryCatch({
    get_earliest_principal_node(cds_sub, marker_gene = "Ly6c2")
  }, error = function(e) {
    message(paste("Warning: Ly6c2 expression too low in", grp, "to find root. Using an arbitrary root instead."))
    # 如果找不到，强制选取图中的第一个主节点作为起点，避免脚本中断
    return(igraph::V(principal_graph(cds_sub)[["UMAP"]])$name[1])
  })
  
  cds_sub <- order_cells(cds_sub, root_pr_nodes = root_node_sub)
  
  # 4.6 可视化当前组
  # 绘制拟时序图
  plot_sub_time <- do.call(plot_cells, c(common_args, list(cds = cds_sub, color_cells_by = "pseudotime"))) + 
    ggtitle(paste("Aorta Macrophages (", grp, "): Pseudotime", sep=""))
  
  # 绘制聚类图
  plot_sub_cluster <- do.call(plot_cells, c(common_args, list(cds = cds_sub, color_cells_by = "cell_annotation"))) + 
    ggtitle(paste("Aorta Macrophages (", grp, "): Clusters", sep=""))
  
  # 4.7 保存图片
  # 替换组名中的特殊字符(如空格)以防文件路径报错
  safe_grp_name <- gsub(" |/|\\\\", "_", grp) 
  
  file_name_time <- paste0("pictures/monocle3_aorta_macrophages_", safe_grp_name, "_pseudotime.png")
  file_name_cluster <- paste0("pictures/monocle3_aorta_macrophages_", safe_grp_name, "_clusters.png")
  
  ggsave(file_name_time, plot_sub_time, width = 6, height = 5)
  ggsave(file_name_cluster, plot_sub_cluster, width = 8, height = 5)
  
  print(paste("Saved:", file_name_time))
}

print("Group-wise trajectory analysis COMPLETED successfully!")