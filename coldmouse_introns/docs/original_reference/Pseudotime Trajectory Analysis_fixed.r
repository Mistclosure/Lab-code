# ==============================================================================
# Script 06 (修复版：解决轨迹断裂，单核-巨噬迁移轨迹)
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

# 2. 提取与标注（核心优化：剔除DC，仅保留目标谱系，消除噪音）
print("Subsetting and annotating target cells (Monocyte/Macrophage only)...")
# 仅保留PBMC中的单核细胞（new_clusters 12），剔除DC（3），消除谱系噪音
pbmc_mono <- subset(pbmc, new_clusters == '12')
pbmc_mono$cell_type <- "Monocyte"
pbmc_mono$Source <- "PBMC_Monocyte"

# 仅保留主动脉中的Macrophages，剔除DC，聚焦核心分化谱系
aorta_mac <- subset(aorta, SingleR.labels == "Macrophages")
aorta_mac$cell_type <- "Macrophages"
aorta_mac$Source <- "Aorta_Macrophages"

# 3. 【核心修复】Seurat锚点整合（彻底解决跨组织批次差异）
print("Seurat Anchor Integration for cross-tissue correction...")
# 构建对象列表，标准化+高可变基因筛选
obj_list <- list(PBMC = pbmc_mono, Aorta = aorta_mac)
obj_list <- lapply(obj_list, function(x){
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# 找跨数据集锚点，整合数据
anchors <- FindIntegrationAnchors(object.list = obj_list, dims = 1:30, k.filter = 50)
merged_obj <- IntegrateData(anchorset = anchors, dims = 1:30)

# 保留原始counts层（Monocle3必须依赖原始counts构建对象）
DefaultAssay(merged_obj) <- "RNA"
merged_obj <- JoinLayers(merged_obj)

# 验证细胞分布
print("=== 合并后细胞分布验证 ===")
print(table(merged_obj$Source, merged_obj$cell_type))

# 4. 构建 Monocle 3 对象
print("Converting to Monocle 3 CDS object...")
expression_matrix <- LayerData(merged_obj, layer = "counts")
cell_metadata <- merged_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                            row.names = rownames(expression_matrix))

cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# 计算Size Factors
cds <- estimate_size_factors(cds)

# ---------------------------------------------------------
# 🌟 修复版Monocle 3 核心流程（解决断裂）
# ---------------------------------------------------------

# 步骤A: PCA预处理（优化维度）
print("Preprocessing with PCA...")
cds <- preprocess_cds(cds, num_dim = 50, method = "PCA")

# 步骤B: 【可选二次对齐】针对残留的微小批次差异
print("Secondary alignment for residual batch effect...")
cds <- align_cds(cds, num_dim = 50, alignment_group = "Source")

# 步骤C: 【核心修复】轨迹优化版UMAP降维（彻底解决簇拆分）
print("Optimized UMAP for continuous trajectory...")
cds <- reduce_dimension(
  cds, 
  reduction_method = "UMAP", 
  preprocess_method = "Aligned", # 严格使用对齐后的结果
  umap.n_neighbors = 50,         # 增大近邻数，捕捉全局谱系连续关系
  umap.min_dist = 0.01,          # 极小化最小距离，让细胞连续聚拢，消除空隙
  umap.metric = "cosine",        # 余弦距离更适配单细胞转录组稀疏数据    
)

# 步骤D: 优化聚类（避免极端分辨率干扰轨迹）
print("Optimized cell clustering...")
cds <- cluster_cells(
  cds, 
  reduction_method = "UMAP", 
  resolution = 0.01,  # 修正极端的1e-3分辨率，适配谱系分析
  k = 50               # 和UMAP近邻数匹配，保证近邻图一致性
)

# 步骤E: 【核心修复】优化轨迹图学习（消除硬拉线）
print("Learning continuous Trajectory Graph...")
cds <- learn_graph(
  cds, 
  use_partition = FALSE,        # 强制全局建图，忽略分区，避免轨迹被截断
  close_loop = FALSE,            # 关闭闭环，符合单核→巨噬的线性分化逻辑
  learn_graph_control = list(
    minimal_branch_len = 3,      # 缩小最小分支长度，保留过渡态细胞的连接
    euclidean_distance_ratio = 10, # 扩大跨簇连接的允许范围
    geodesic_distance_ratio = 1/5  # 放宽近邻连接的限制
  )
)

# 步骤F: 基于Ly6c2高表达单核细胞设定拟时序起点（完全保留你的需求）
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
print("Generating optimized Plots...")

# 来源分布图（看PBMC和Aorta是否连续连接）
plot_source_global <- plot_cells(
  cds, 
  color_cells_by = "Source", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE
) + 
  ggtitle("Fixed Trajectory: PBMC Monocyte to Aorta Macrophage")

# 全局伪时间图（看分化方向是否连续）
plot_pseudotime_global <- plot_cells(
  cds, 
  color_cells_by = "pseudotime", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE
) + 
  ggtitle("Continuous Pseudotime Trajectory")

# 细胞类型分布图（验证谱系连续性）
plot_celltype <- plot_cells(
  cds, 
  color_cells_by = "cell_type", 
  label_cell_groups = FALSE, 
  label_leaves = FALSE, 
  label_branch_points = FALSE
) + 
  ggtitle("Cell Type Distribution")

# 分组对比图（如果有Group分组）
if("Group" %in% colnames(colData(cds))){
  plot_pseudotime_split <- plot_cells(
    cds, 
    color_cells_by = "pseudotime", 
    label_cell_groups = FALSE, 
    label_leaves = FALSE, 
    label_branch_points = FALSE
  ) + 
    facet_wrap(~Group) + 
    ggtitle("Pseudotime Comparison: Cold vs RT")
}

# 6. 保存
print("Saving outputs...")
if(!dir.exists("files")) dir.create("files")
if(!dir.exists("pictures")) dir.create("pictures")

qsave(cds, "files/script_06_monocle3_fixed_trajectory.qs")

ggsave("pictures/monocle3_fixed_source.png", plot_source_global, width = 7, height = 6)
ggsave("pictures/monocle3_fixed_pseudotime.png", plot_pseudotime_global, width = 7, height = 6)
ggsave("pictures/monocle3_fixed_celltype.png", plot_celltype, width = 7, height = 6)
if(exists("plot_pseudotime_split")){
  ggsave("pictures/monocle3_fixed_split.png", plot_pseudotime_split, width = 12, height = 6)
}

print("Script 06 FIXED VERSION ALL COMPLETED successfully!")