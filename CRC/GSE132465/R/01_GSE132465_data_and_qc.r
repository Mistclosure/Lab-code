# ==============================================================================
# GSE132465 专属版本：单一大矩阵读取、全局质控去核糖体、按Patient分批去双胞
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(patchwork)
library(data.table)

if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/')

# ------------------------------------------------------------------------------
# 1. 庞大 UMI 矩阵数据读取与 Seurat 对象创建
# ------------------------------------------------------------------------------
print("🚀 步骤1/4: 开始读取 10X 矩阵与注释数据...")

# 读取矩阵
counts_matrix <- fread("GSE132465_10X_UMI_matrix.txt", header = TRUE, sep = "\t", check.names = FALSE)
counts_matrix <- as.data.frame(counts_matrix)
rownames(counts_matrix) <- counts_matrix[[1]] # 第一列设为行名
counts_matrix <- counts_matrix[,-1]          # 删掉原来的第一列

# 读取注释
cell_annotation <- fread("GSE132465_10X_cell_annotation.txt", header = TRUE, sep = "\t")
cell_annotation <- as.data.frame(cell_annotation)
rownames(cell_annotation) <- cell_annotation[[1]]
cell_annotation <- cell_annotation[,-1]

# 创建 Seurat 对象
sc_combined <- CreateSeuratObject(counts = counts_matrix, meta.data = cell_annotation)

# 【优化 1：极限节省内存】一旦数据挂载进 Seurat 对象，立刻清除 R 环境里庞大的明文矩阵变量！
rm(counts_matrix, cell_annotation)
gc()

# 计算线粒体比例
if (any(grepl("^MT-", rownames(sc_combined)))) {
  sc_combined[["percent.mt"]] <- PercentageFeatureSet(sc_combined, pattern = "^MT-")
} else {
  sc_combined[["percent.mt"]] <- PercentageFeatureSet(sc_combined, pattern = "^mt-")
}

print(paste("   ✅ 数据加载完成，初始细胞总数:", ncol(sc_combined)))

# ------------------------------------------------------------------------------
# 2. 综合质控：基础指标 + 核糖体彻底处理
# ------------------------------------------------------------------------------
print("🚀 步骤2/4: 正在进行综合质控 (过滤指标并彻底移除核糖体基因)...")

# 基础过滤
sc_combined <- subset(sc_combined, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 15)

# 寻找核糖体基因并移除
sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = "^[Rr][Pp][SsLl]")
ribo_genes <- grep("^[Rr][Pp][SsLl]", rownames(sc_combined), value = TRUE)
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]

# 严格保留非核糖体基因
sc_combined <- subset(sc_combined, features = non_ribo_genes)

print(paste("   ✅ 过滤与移除核糖体基因完成，当前基因数:", nrow(sc_combined), "| 细胞数:", ncol(sc_combined)))

# ------------------------------------------------------------------------------
# 3. 按 Patient 拆分并独立去双胞 (scDblFinder)
# ------------------------------------------------------------------------------
print("🚀 步骤3/4: 按 metadata 中的 Patient 拆分并执行 scDblFinder 去双胞...")

# 【防呆检查】：确保你的 metadata 里确实有一列叫 Patient，否则提前报错
if (!"Patient" %in% colnames(sc_combined@meta.data)) {
  stop("❌ 错误：在 metadata 中未找到 'Patient' 这一列！请检查 GSE132465_10X_cell_annotation.txt 的列名，或手动重命名列。")
}

# 按 Patient 进行拆分
sc_split <- SplitObject(sc_combined, split.by = "Patient")

# 【优化 2：内存释放】拆分完 list 后，抹去包含全部细胞的大合集，防止双倍内存占用
rm(sc_combined)
gc()

run_standard_pipeline <- function(obj, patient_name) {
  print(paste(">>> 正在处理患者批次:", patient_name, "| 质控后细胞数:", ncol(obj)))
  if(ncol(obj) < 50) { return(NULL) }
  
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5) 
  
  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    # 独立计算每个患者的双细胞
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class
    
    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    print(paste("   [scDblFinder] 患者", patient_name, "发现双细胞:", n_dbl, "个 (占比", round(n_dbl/ncol(obj)*100, 2), "%)"))
    
    # 仅保留单细胞
    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    
  }, error = function(e) {
    message(paste("   ⚠️ scDblFinder 运行警告:", e$message, "- 患者", patient_name, "将保留所有细胞继续"))
  })
  
  return(obj)
}

# 循环运行，此处以 sapply/lapply 处理所有批次
sc_clean_list <- lapply(names(sc_split), function(x) {
  res <- run_standard_pipeline(sc_split[[x]], x)
  # 【优化 3】每次循环后回收一次内存
  gc()
  return(res)
})
sc_clean_list <- sc_clean_list[!sapply(sc_clean_list, is.null)]

# ------------------------------------------------------------------------------
# 4. 重新合并与保存
# ------------------------------------------------------------------------------
print("🚀 步骤4/4: 合并纯净数据并保存...")

sc_final <- merge(sc_clean_list[[1]], y = sc_clean_list[2:length(sc_clean_list)])
sc_final <- JoinLayers(sc_final)

qsave(sc_final, "sc_combined_qc_Cleaned.qs")

print(paste("🎉 全部运行完毕！最终保留高质量单细胞总数:", ncol(sc_final)))
print("📊 各患者(Patient)的最终保留细胞数统计：")
print(table(sc_final$Patient))

# 结束清理环境
rm(sc_clean_list)
gc()