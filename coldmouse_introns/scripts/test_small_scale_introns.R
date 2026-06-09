# ==============================================================================
# test_small_scale_introns.R: 小规模测试脚本
# ==============================================================================
# 只读取 A1.loom，验证 loom 结构、Seurat 对象、QC 分布、小规模聚类、monocle2 转换。
# 输出保存到 results/test/，不覆盖完整结果。
# ==============================================================================

# 加载工具
PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))
source(file.path(PROJECT_HOME, "R", "utils_qc.R"))
source(file.path(PROJECT_HOME, "R", "utils_seurat.R"))
source(file.path(PROJECT_HOME, "R", "utils_monocle2.R"))

config <- load_config()

# 加载 R 包
library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(hdf5r)
library(patchwork)
library(scDblFinder)
library(SingleCellExperiment)

# 创建测试输出目录
test_dir <- get_results_dir(config, "test")

# ------------------------------------------------------------------------------
# 1. 检查 loom 文件结构
# ------------------------------------------------------------------------------
print("=== 测试1: 检查 loom 文件结构 ===")
loom_file <- file.path(config$paths$rawdata_dir, "A1.loom")
loom_info <- check_loom_structure(loom_file)

cat("Loom 文件:", loom_file, "\n")
cat("细胞数:", loom_info$n_cells, "\n")
cat("基因数:", loom_info$n_genes, "\n")
cat("矩阵维度:", loom_info$matrix_dims, "\n")
cat("是否 gene x cell:", loom_info$is_gene_x_cell, "\n")
cat("前10个基因:", head(loom_info$genes, 10), "\n")
cat("前10个barcode:", head(loom_info$barcodes, 10), "\n")

# 保存 loom 结构信息
loom_summary <- data.frame(
  file = loom_file,
  n_cells = loom_info$n_cells,
  n_genes = loom_info$n_genes,
  is_gene_x_cell = loom_info$is_gene_x_cell
)
write.csv(loom_summary, file.path(test_dir, "test_loom_structure.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 2. 读取 A1.loom 并创建 Seurat 对象
# ------------------------------------------------------------------------------
print("=== 测试2: 读取 A1.loom 并创建 Seurat 对象 ===")
counts <- read_introns_loom(loom_file, sample = "A1", chunk_cells = 2000)

cat("Counts 矩阵维度:", dim(counts), "\n")
cat("非零元素:", nnzero(counts), "\n")
cat("稀疏度:", round(1 - nnzero(counts) / (nrow(counts) * ncol(counts)), 4), "\n")

# 创建 Seurat 对象
sc_obj <- CreateSeuratObject(counts = counts, project = "A1", min.cells = 3, min.features = 200)
cat("Seurat 对象细胞数:", ncol(sc_obj), "\n")
cat("Seurat 对象基因数:", nrow(sc_obj), "\n")
cat("Assay 名称:", names(sc_obj@assays), "\n")
cat("Metadata 列:", colnames(sc_obj@meta.data), "\n")

# 添加 metadata
sc_obj$Tissue <- "Aorta"
sc_obj$Group <- "RT_25C"
sc_obj$Condition <- "NonCold"
sc_obj <- calculate_qc_metrics(sc_obj)

cat("添加 metadata 后列:", colnames(sc_obj@meta.data), "\n")
cat("nCount_RNA 范围:", range(sc_obj$nCount_RNA), "\n")
cat("nFeature_RNA 范围:", range(sc_obj$nFeature_RNA), "\n")
cat("percent.mt 范围:", range(sc_obj$percent.mt), "\n")
cat("percent.rb 范围:", range(sc_obj$percent.rb), "\n")

# ------------------------------------------------------------------------------
# 3. QC 分布检查
# ------------------------------------------------------------------------------
print("=== 测试3: QC 分布检查 ===")
qc_quantiles <- calculate_qc_quantiles(sc_obj)
print(qc_quantiles)
write.csv(qc_quantiles, file.path(test_dir, "test_qc_quantiles_A1.csv"), row.names = FALSE)

# QC 小提琴图
p_vln <- VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"),
                  ncol = 4, pt.size = 0)
ggsave(file.path(test_dir, "test_qc_violin_A1.png"), plot = p_vln, width = 12, height = 6, dpi = 150)

# QC 散点图
p_scatter1 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p_scatter2 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
p_scatter <- p_scatter1 + p_scatter2
ggsave(file.path(test_dir, "test_qc_scatter_A1.png"), plot = p_scatter, width = 12, height = 5, dpi = 150)

# ------------------------------------------------------------------------------
# 4. 小规模 NormalizeData, FindVariableFeatures, ScaleData, PCA, UMAP
# ------------------------------------------------------------------------------
print("=== 测试4: 小规模 Seurat 流程 ===")
sc_obj <- JoinLayers(sc_obj)
sc_obj <- NormalizeData(sc_obj)
sc_obj <- FindVariableFeatures(sc_obj, selection.method = "vst", nfeatures = 2000)
sc_obj <- ScaleData(sc_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
sc_obj <- RunPCA(sc_obj, npcs = 30, verbose = FALSE)
sc_obj <- FindNeighbors(sc_obj, dims = 1:20)
sc_obj <- FindClusters(sc_obj, resolution = 0.5)
sc_obj <- RunUMAP(sc_obj, dims = 1:20)

cat("PCA 维度:", dim(Embeddings(sc_obj, "pca")), "\n")
cat("UMAP 维度:", dim(Embeddings(sc_obj, "umap")), "\n")
cat("聚类数:", length(unique(Idents(sc_obj))), "\n")
cat("聚类分布:\n")
print(table(Idents(sc_obj)))

# UMAP 图
p_umap <- DimPlot(sc_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("A1 Test UMAP")
ggsave(file.path(test_dir, "test_umap_A1.png"), plot = p_umap, width = 8, height = 6, dpi = 150)

# ------------------------------------------------------------------------------
# 5. 测试 Seurat v5 -> monocle2 转换
# ------------------------------------------------------------------------------
print("=== 测试5: Seurat -> Monocle2 转换 ===")
tryCatch({
  library(monocle)
  cds <- seurat_to_monocle2(sc_obj)
  cat("Monocle2 CellDataSet 细胞数:", ncol(cds), "\n")
  cat("Monocle2 CellDataSet 基因数:", nrow(cds), "\n")
  cat("Phenodata 列:", colnames(pData(cds)), "\n")
  cat("Featuredata 列:", colnames(fData(cds)), "\n")
  saveRDS(cds, file.path(test_dir, "test_monocle2_cds_A1.rds"))
  print("Monocle2 转换成功!")
}, error = function(e) {
  message("Monocle2 转换失败: ", e$message)
})

# 保存测试 Seurat 对象
saveRDS(sc_obj, file.path(test_dir, "test_seurat_A1.rds"))

print("=== 所有测试完成 ===")
print(paste("测试输出保存在:", test_dir))
