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
cds_mono <- reduce_dimension(cds_mono, reduction_method = "UMAP")
cds_mono <- cluster_cells(cds_mono, reduction_method = "UMAP", resolution = 1e-3)

# 4. 学习轨迹图
print("Learning Trajectory Graph...")
cds_mono <- learn_graph(cds_mono, use_partition = FALSE, 
                        learn_graph_control = list(minimal_branch_len = 5))

# 5. 设定拟时序起点并排序
print("Ordering cells...")
# 同样假设 RT_25C 为起点
root_candidates <- colnames(cds_mono)[colData(cds_mono)$Group == "RT_25C"]
if(length(root_candidates) == 0) root_candidates <- colnames(cds_mono)
cds_mono <- order_cells(cds_mono, root_cells = root_candidates)

# 6. 可视化
print("Generating Plots...")
plot_group <- plot_cells(cds_mono, color_cells_by = "Group", 
                         label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("PBMC Monocytes: Trajectory by Group")

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

qsave(cds_mono, "files/script_06_03_monocle3_pbmc_monocytes.qs")
ggsave("pictures/monocle3_06_03_group.pdf", plot_group, width = 6, height = 5)
ggsave("pictures/monocle3_06_03_pseudotime.pdf", plot_time, width = 6, height = 5)
ggsave("pictures/monocle3_06_03_split.pdf", plot_split, width = 10, height = 5)

print("Script 06.03 COMPLETED successfully!")