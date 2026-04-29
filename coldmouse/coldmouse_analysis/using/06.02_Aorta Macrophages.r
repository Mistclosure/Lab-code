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

# 3. 降维与聚类 (不再需要跨组织对齐)
print("Preprocessing and dimension reduction...")
cds_mac <- preprocess_cds(cds_mac, num_dim = 30)
cds_mac <- reduce_dimension(cds_mac, reduction_method = "UMAP")
cds_mac <- cluster_cells(cds_mac, reduction_method = "UMAP", resolution = 1e-3)

# 4. 学习轨迹图
print("Learning Trajectory Graph...")
cds_mac <- learn_graph(cds_mac, use_partition = FALSE, 
                       learn_graph_control = list(minimal_branch_len = 5))

# 5. 设定拟时序起点并排序
print("Ordering cells...")
# 假设 RT_25C 是基线/未激活状态，将其设为起点。如果找不到 RT_25C，则回退到随机起点以防报错
root_candidates <- colnames(cds_mac)[colData(cds_mac)$Group == "RT_25C"]
if(length(root_candidates) == 0) root_candidates <- colnames(cds_mac) 
cds_mac <- order_cells(cds_mac, root_cells = root_candidates)

# 6. 可视化
print("Generating Plots...")
plot_group <- plot_cells(cds_mac, color_cells_by = "Group", 
                         label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("Aorta Macrophages: Trajectory by Group")

plot_time <- plot_cells(cds_mac, color_cells_by = "pseudotime", 
                        label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("Aorta Macrophages: Pseudotime")

plot_split <- plot_cells(cds_mac, color_cells_by = "pseudotime", 
                         label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  facet_wrap(~Group) + 
  ggtitle("Aorta Macrophages: Pseudotime Split by Group")

# 7. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

qsave(cds_mac, "files/script_06_02_monocle3_aorta_macrophages.qs")
ggsave("pictures/monocle3_06_02_group.pdf", plot_group, width = 6, height = 5)
ggsave("pictures/monocle3_06_02_pseudotime.pdf", plot_time, width = 6, height = 5)
ggsave("pictures/monocle3_06_02_split.pdf", plot_split, width = 10, height = 5)

print("Script 06.02 COMPLETED successfully!")