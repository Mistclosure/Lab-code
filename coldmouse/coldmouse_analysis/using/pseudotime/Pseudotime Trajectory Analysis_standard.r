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
library(harmony)

# 1. 读取数据
print("Loading data...")
pbmc <- qread('pbmc_corrected.qs')
aorta <- qread('aorta_corrected.qs')

# 2. 提取与标注 (修改版：分别定义PBMC中的Monocyte和DC)
print("Subsetting and annotating PBMC cells...")
#先提取PBMC中new_clusters为3和12的细胞
pbmc_sub <- subset(pbmc, new_clusters %in% c('3', '12'))
#定义细胞类型：12→Monocyte，3→Dendritic cells
pbmc_sub$cell_type <- ifelse(pbmc_sub$new_clusters == '12', 
                              "Monocyte", "Dendritic cells")
# 定义Source来源
#pbmc_sub$Source <- paste0("PBMC_", pbmc_sub$celltype)
pbmc_sub = subset(pbmc_sub, new_clusters == '12')
print("Subsetting Aorta cells...")
# 提取主动脉中的Macrophages和Dendritic cells
aorta_target <- subset(aorta, seurat_clusters %in% c('1','6','10','14'))
# 统一细胞类型列名（方便后续合并分析）
aorta_target$cell_type <- aorta_target$SingleR.labels
# 定义Source来源
#aorta_target$Source <- paste0("Aorta_", aorta_target$celltype)

# 3. 合并与预处理
print("Merging objects...")
# 合并两个对象，添加细胞前缀避免重复
merged_obj <- merge(pbmc_sub, y = aorta_target, 
                    add.cell.ids = c("PBMC", "Aorta"))
merged_obj = subset(merged_obj,Group == 'TN_30C',invert = TRUE)
merged_obj$celllabels <- sub(".*:", "", merged_obj$celltype)
merged_obj$Source <- merged_obj$celllabels
# 合并图层（Seurat v5必需步骤）
merged_obj <- JoinLayers(merged_obj)

# 可选：验证合并后的细胞类型分布
print(table(merged_obj$Source, merged_obj$cell_type))

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
cds <- align_cds(cds, num_dim = 30, alignment_group = "Source",alignment_method = "harmony")

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

# 步骤 F: 设定拟时序起点并排序（修改为基于Ly6c2高表达单核细胞）
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

root_node <- get_earliest_principal_node(cds, marker_gene = "Ly6c2")
cds <- order_cells(cds, root_pr_nodes = root_node)

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

#qsave(cds, "files/script_06_monocle3_native_standard.qs")

ggsave("pictures/monocle3_standard_source01.png", plot_source_global, width = 8, height = 5)
ggsave("pictures/monocle3_standard_pseudotime01.png", plot_pseudotime_global, width = 6, height = 5)
ggsave("pictures/monocle3_standard_split01.png", plot_pseudotime_split, width = 10, height = 5)

print("Script 06 ALL COMPLETED successfully!")