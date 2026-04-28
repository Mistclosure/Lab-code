# ==============================================================================
# Script 06 (Ultra-Robust Version): Seurat v5 & Monocle 3 Trajectory
# Purpose: Research Monocyte Migration (PBMC -> Aorta Monocytes -> Macrophages)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse/')
library(qs)
library(Seurat)
library(monocle3)
library(ggplot2)
library(dplyr)

# 1. 读取数据
print("Loading data...")
pbmc <- qread('pbmc_corrected.qs')
sc_by_tissue <- qread('sc_by_tissue_louvain.qs')
aorta <- sc_by_tissue[['Aorta']]

# 2. 提取与标注
print("Subsetting Monocytes and Macrophages...")
pbmc_mono <- subset(pbmc, new_clusters == "12") 
pbmc_mono$Source <- "PBMC_Monocyte"

aorta_target <- subset(aorta, SingleR.labels %in% c("Macrophages"))
aorta_target$Source <- paste0("Aorta_", aorta_target$SingleR.labels)

# 3. 合并与经典预处理
print("Merging and Preprocessing...")
merged_obj <- merge(pbmc_mono, y = aorta_target, add.cell.ids = c("PBMC", "Aorta"))
merged_obj <- JoinLayers(merged_obj)

merged_obj <- NormalizeData(merged_obj) %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = FALSE)

# 4. 构建全局 Monocle 3 对象
print("Converting to Monocle 3...")
expression_matrix <- LayerData(merged_obj, layer = "counts")
cell_metadata <- merged_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                            row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)
cds <- estimate_size_factors(cds)

print("Injecting PCA reliably...")
reducedDims(cds)[["PCA"]] <- Embeddings(merged_obj, reduction = "pca")
pca_variance <- merged_obj[["pca"]]@stdev^2
attr(reducedDims(cds)[["PCA"]], "proportion_var_expl") <- pca_variance / sum(pca_variance)

# ---------------------------------------------------------
# 🌟 关键修复点：使用 align_cds 抹平组织间的巨大转录鸿沟
# ---------------------------------------------------------
print("Aligning PBMC and Aorta (Stitching the islands)...")
# 这步会基于我们注入的 PCA，运用 MNN 算法将 PBMC 和 Aorta 的边缘强行拼接
cds <- align_cds(cds, num_dim = 30, alignment_group = "Source")

# 5. 基于【对齐后】的坐标学习全局轨迹
print("Learning Global Trajectory Graph...")
# 注意：preprocess_method 改为了 "Aligned"
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "Aligned")

# 顺便优化一下聚类参数，防止左侧（拟时序起点）出现过度密集的乱麻状连线
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3) 
cds <- learn_graph(cds, use_partition = FALSE, 
                   learn_graph_control = list(minimal_branch_len = 10)) # 砍掉太细碎的分支

# 6. 设定拟时序起点 (强化版)
print("Ordering cells based on PBMC Monocytes...")
# ⭐️ 核心修正：直接从 cds 对象中提取条形码，绝对防止名称不匹配
monocyte_cells <- colnames(cds)[colData(cds)$Source == "PBMC_Monocyte"]
cds <- order_cells(cds, root_cells = monocyte_cells)

# 7. 可视化
print("Generating Plots...")
plot_source_global <- plot_cells(cds, color_cells_by = "Source", 
                                 label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("Global Trajectory: PBMC to Aorta (Mono & Mac)")

plot_pseudotime_global <- plot_cells(cds, color_cells_by = "pseudotime", 
                                     label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  ggtitle("Global Pseudotime Trajectory")

plot_pseudotime_split <- plot_cells(cds, color_cells_by = "pseudotime", 
                                    label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE) + 
  facet_wrap(~Group) + 
  ggtitle("Pseudotime Comparison: Cold vs RT")

# 8. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
qsave(cds, "files/script_06_monocle3_robust_injected.qs")

ggsave("files/monocle3_global_source.pdf", plot_source_global, width = 6, height = 5)
ggsave("files/monocle3_global_pseudotime.pdf", plot_pseudotime_global, width = 6, height = 5)
ggsave("files/monocle3_split_by_group.pdf", plot_pseudotime_split, width = 10, height = 5)

print("Script 06 ALL COMPLETED successfully!")
