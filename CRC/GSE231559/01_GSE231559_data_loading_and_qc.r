# ==============================================================================
# 综合版脚本：多样本自动读取、元数据注入、核糖体移除与自动化去双胞
# ==============================================================================

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(patchwork)

if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

# ------------------------------------------------------------------------------
# 1. 设置路径与获取样本列表
# ------------------------------------------------------------------------------
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE231559')
data_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE231559/rawdata" 
sample_ids <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)
sample_ids <- sample_ids[grepl("^[A-Z][0-9]+[A-Z]", sample_ids)]

# ------------------------------------------------------------------------------
# 2. 数据读取与 Seurat 对象创建 (包含 Patient 和 Type 提取)
# ------------------------------------------------------------------------------
sc_list <- list()

print("🚀 步骤1/4: 开始读取 10X 矩阵数据并注入元数据...")

for (sample in sample_ids) {
  folder_path <- file.path(data_dir, sample)
  
  tryCatch({
    counts <- Read10X(data.dir = folder_path)
    # 使用文件夹名称作为 project，这会自动存入 metadata 的 orig.ident
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    if (ncol(sc_obj) == 0) next 
    
    sc_obj$Patient <- str_extract(sample, "^[A-Z][0-9]+")
    
    if (grepl("T", sample)) {
      sc_obj$Type <- "Tumor"
    } else if (grepl("N", sample)) {
      sc_obj$Type <- "Normal"
    } else {
      sc_obj$Type <- "Unknown"
    }
    
    if (any(grepl("^MT-", rownames(sc_obj)))) {
      sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
    } else {
      sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
    }
    
    sc_list[[sample]] <- sc_obj
    print(paste("   ✅ 样本加载完成:", sample, "| 患者:", unique(sc_obj$Patient), "| 类型:", unique(sc_obj$Type)))
    
  }, error = function(e) { 
    message(paste("   ❌ 处理样本", sample, "时出错:", e$message)) 
  })
}

sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = names(sc_list))
sc_combined <- JoinLayers(sc_combined) 

# ------------------------------------------------------------------------------
# 3. 综合质控：基础指标 + 核糖体处理
# ------------------------------------------------------------------------------
print("🚀 步骤2/4: 正在进行综合质控 (过滤指标并彻底移除核糖体基因)...")

sc_combined <- subset(sc_combined, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 15)

sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = "^[Rr][Pp][SsLl]")
ribo_genes <- grep("^[Rr][Pp][SsLl]", rownames(sc_combined), value = TRUE)
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]

sc_combined <- subset(sc_combined, features = non_ribo_genes)
print(paste("   ✅ 移除核糖体基因完成，当前基因数:", nrow(sc_combined)))

# ------------------------------------------------------------------------------
# 4. 按样本 (orig.ident) 拆分并独立去双胞 (scDblFinder)
# ------------------------------------------------------------------------------
print("🚀 步骤3/4: 按样本 (Sample/orig.ident) 拆分并执行 scDblFinder 去双胞...")

# --- 核心改动：按 orig.ident 拆分，确保每个样本独立去双胞 ---
sc_split <- SplitObject(sc_combined, split.by = "orig.ident")

run_standard_pipeline <- function(obj, sample_name) {
  print(paste(">>> 正在处理样本:", sample_name, "| 质控后细胞数:", ncol(obj)))
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
    # 独立计算每个样本的双细胞
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class
    
    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    print(paste("   [scDblFinder] 样本", sample_name, "发现双细胞:", n_dbl, "个 (占比", round(n_dbl/ncol(obj)*100, 2), "%)"))
    
    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    
  }, error = function(e) {
    message(paste("   ⚠️ scDblFinder 运行警告:", e$message, "- 样本", sample_name, "将保留所有细胞继续"))
  })
  
  return(obj)
}

sc_clean_list <- lapply(names(sc_split), function(x) {
  run_standard_pipeline(sc_split[[x]], x)
})
sc_clean_list <- sc_clean_list[!sapply(sc_clean_list, is.null)]

# ------------------------------------------------------------------------------
# 5. 重新合并与保存
# ------------------------------------------------------------------------------
print("🚀 步骤4/4: 合并纯净数据并保存...")

sc_final <- merge(sc_clean_list[[1]], y = sc_clean_list[2:length(sc_clean_list)])
sc_final <- JoinLayers(sc_final)

qsave(sc_final, "sc_combined_qc_Cleaned.qs")

print(paste("🎉 全部运行完毕！最终保留高质量单细胞总数:", ncol(sc_final)))
print("📊 样本及患者分布统计：")
print(table(sc_final$Patient, sc_final$Type))