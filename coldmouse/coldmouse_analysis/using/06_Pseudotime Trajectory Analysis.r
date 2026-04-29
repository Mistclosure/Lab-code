# ==============================================================================
# Script 06 (Standard Monocle 3 Native Workflow)
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

# 3. 合并与经典预处理 (Seurat阶段仅做合并)
print("Merging objects...")
merged_obj <- merge(pbmc_mono, y = aorta_target, add.cell.ids = c("PBMC", "Aorta"))
merged_obj <- JoinLayers(merged_obj)

# 4. 构建全局 Monocle 3 对象
print("Converting to Monocle 3...")
expression_matrix <- LayerData(merged_obj, layer = "counts")
cell_metadata <- merged_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                            row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# 计算 Size Factors (Monocle 的归一化步骤)
cds <- estimate_size_factors(cds)

# ---------------------------------------------------------
# 🌟 Monocle 3 原生核心分析流程开始
# ---------------------------------------------------------

# 步骤 A: 原生 PCA 预处理
print("Preprocessing natively in Monocle 3 (PCA)...")
cds <- preprocess_cds(cds, num_dim = 30)

# 步骤 B: 批次效应/组织差异对齐 (使用 MNN)
print("Aligning PBMC and Aorta...")
cds <- align_cds(cds, num_dim = 30, alignment_group = "Source")

# 步骤 C: UMAP 降维 (基于对齐后的结果)
# 注意：这里必须使用 UMAP，后续的轨迹推断强依赖于此
print("Reducing dimensions (UMAP)...")
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = "Aligned")

# 步骤 D: 聚类
print("Clustering cells...")
cds <- cluster_cells(cds, reduction_method = "UMAP", resolution = 1e-3) 

# 步骤 E: 学习轨迹图
print("Learning Trajectory Graph...")
cds <- learn_graph(cds, use_partition = FALSE, 
                   learn_graph_control = list(minimal_branch_len = 10))

# 步骤 F: 设定拟时序起点并排序细胞
print("Ordering cells based on PBMC Monocytes...")
# 提取 PBMC_Monocyte 的细胞条形码作为拟时序起点 (Root)
monocyte_cells <- colnames(cds)[colData(cds)$Source == "PBMC_Monocyte"]
cds <- order_cells(cds, root_cells = monocyte_cells)

# ---------------------------------------------------------
# 5. 可视化
# ---------------------------------------------------------
print("Generating Plots...")

# 来源分布图 (看 PBMC 和 Aorta 是否连上了)
plot_source_global <- plot_cells(cds, color_cells_by = "Source", 
                                 label_cell_groups = FALSE, 
                                 label_leaves = FALSE, 
                                 label_branch_points = FALSE) + 
  ggtitle("Global Trajectory: PBMC to Aorta (Mono & Mac)")

# 全局伪时间图 (看演化方向)
plot_pseudotime_global <- plot_cells(cds, color_cells_by = "pseudotime", 
                                     label_cell_groups = FALSE, 
                                     label_leaves = FALSE, 
                                     label_branch_points = FALSE) + 
  ggtitle("Global Pseudotime Trajectory")

# 如果你的 metadata 中有 Group (例如 Cold 和 RT)，可以用这个图拆分对比
plot_pseudotime_split <- plot_cells(cds, color_cells_by = "pseudotime", 
                                    label_cell_groups = FALSE, 
                                    label_leaves = FALSE, 
                                    label_branch_points = FALSE) + 
  facet_wrap(~Group) + 
  ggtitle("Pseudotime Comparison: Cold vs RT")

# 6. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

qsave(cds, "files/script_06_monocle3_native_standard.qs")

ggsave("pictures/monocle3_standard_source.pdf", plot_source_global, width = 6, height = 5)
ggsave("pictures/monocle3_standard_pseudotime.pdf", plot_pseudotime_global, width = 6, height = 5)
ggsave("pictures/monocle3_standard_split.pdf", plot_pseudotime_split, width = 10, height = 5)

print("Script 06 ALL COMPLETED successfully!")