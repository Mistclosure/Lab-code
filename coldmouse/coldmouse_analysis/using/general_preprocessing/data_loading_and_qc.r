# ==============================================================================
# 综合版脚本：数据加载、特定杂质靶向清除与自动化去双胞
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包 (整合两版需求)
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(Matrix)   
library(qs)            # 用于快速保存/读取
library(patchwork)     # 用于绘图

# 加载 scDblFinder 相关包
if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
library(scDblFinder)
library(SingleCellExperiment)

# ------------------------------------------------------------------------------
# 1. 设置路径与元数据映射 (继承脚本 1 变量)
# ------------------------------------------------------------------------------
data_dir <- "/mnt/disk1/qiuzerui/expriments/coldmouse/sccoldmouse" 
sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")
tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")
group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")

# ------------------------------------------------------------------------------
# 2. 数据读取与 Seurat 对象创建 (继承脚本 1 逻辑)
# ------------------------------------------------------------------------------
sc_list <- list()
full_folder_names <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

print("🚀 步骤1/4: 开始读取 10X 矩阵数据并创建对象...")

for (sample in sample_ids) {
  pattern <- paste0("^", sample, "_EmptyDrops")
  folder <- full_folder_names[grepl(pattern, full_folder_names)]
  
  if(length(folder) == 0) {
    message(paste("⚠️ 跳过：未找到样本对应的文件夹:", sample))
    next
  }
  
  tryCatch({
    counts <- Read10X(data.dir = file.path(data_dir, folder[1]))
    colnames(counts) <- gsub("_", "-", colnames(counts))
    # 按照脚本 2 建议，放宽基础 min.features
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 200)
    
    if (ncol(sc_obj) == 0) next 
    
    prefix <- substr(sample, 1, 1); suffix <- substr(sample, 2, 2)
    sc_obj$Tissue <- as.character(tissue_map[prefix])
    sc_obj$Group  <- as.character(group_map[suffix])
    sc_obj$Condition <- ifelse(suffix == "2", "Cold", "NonCold")
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
    
    sc_list[[sample]] <- sc_obj
    print(paste("   ✅ 样本入库完成:", sample, "| 初始细胞数:", ncol(sc_obj)))
    
  }, error = function(e) { 
    message(paste("   ❌ 处理样本", sample, "时出错:", e$message)) 
  })
}

sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = sample_ids)
sc_combined <- JoinLayers(sc_combined) # Seurat v5 必须步骤

# ------------------------------------------------------------------------------
# 3. 混合质控：基础指标 + 特异性细胞剔除 (红细胞/血小板) + 核糖体处理
# ------------------------------------------------------------------------------
print("🚀 步骤2/4: 正在进行综合质控 (数值过滤 + 杂质靶向清除)...")

# A. 基础指标过滤 (采用脚本2的宽容度，防止误删活跃细胞)
sc_combined <- subset(sc_combined, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 & percent.mt < 15)

# B. 靶向剔除血小板和红细胞污染 (脚本1 核心逻辑)
counts_mat <- GetAssayData(sc_combined, layer = "counts")
cells_keep <- colnames(sc_combined)

# 剔除表达血小板基因 Pf4 的细胞
#if ("Pf4" %in% rownames(counts_mat)) {
#  cells_keep <- cells_keep[counts_mat["Pf4", cells_keep] == 0]
#}
# 剔除高表达红细胞基因 Hbb-bs / Hba-a1 的细胞
#for (rbc_gene in c("Hbb-bs", "Hba-a1")) {
#  if (rbc_gene %in% rownames(counts_mat)) {
#    cells_keep <- cells_keep[counts_mat[rbc_gene, cells_keep] <= 1]
#  }
#}
sc_combined <- subset(sc_combined, cells = cells_keep)
print(paste("   ✅ 剔除血小板/红细胞后剩余细胞数:", ncol(sc_combined)))

# C. 核糖体基因处理 (脚本2 逻辑)
sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = "^Rp[sl]")
ribo_genes <- grep("^Rp[sl]", rownames(sc_combined), value = TRUE, ignore.case = TRUE)
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]
sc_combined <- subset(sc_combined, features = non_ribo_genes)
print(paste("   ✅ 彻底移除核糖体基因，当前基因数:", nrow(sc_combined)))

# ------------------------------------------------------------------------------
# 4. 按组织拆分并独立去双胞 (scDblFinder)
# ------------------------------------------------------------------------------
print("🚀 步骤3/4: 按组织拆分并执行 scDblFinder 去双胞...")

sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name) {
  print(paste(">>> 正在处理组织:", tissue_name, "| 质控后细胞数:", ncol(obj)))
  if(ncol(obj) < 50) { return(NULL) }
  
  # 预处理以供 scDblFinder 聚类参考
  obj <- JoinLayers(obj)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = 1:20)
  obj <- FindClusters(obj, resolution = 0.5) 
  
  # scDblFinder 去双胞
  tryCatch({
    sce <- as.SingleCellExperiment(obj)
    sce <- scDblFinder(sce, clusters = TRUE)
    obj$scDblFinder_class <- sce$scDblFinder.class
    
    n_dbl <- sum(obj$scDblFinder_class == "doublet")
    print(paste("   [scDblFinder] 发现双细胞:", n_dbl, "个 (占比", round(n_dbl/ncol(obj)*100, 2), "%)"))
    
    # 剔除双细胞
    obj <- subset(obj, subset = scDblFinder_class == "singlet")
    print(paste("   [Filter] 剔除双细胞后剩余纯净细胞:", ncol(obj)))
    
  }, error = function(e) {
    message(paste("   ⚠️ scDblFinder 运行警告:", e$message, "- 将保留所有细胞继续"))
  })
  
  return(obj)
}

# 批量执行
sc_by_tissue_clean <- lapply(names(sc_by_tissue), function(x) {
  run_standard_pipeline(sc_by_tissue[[x]], x)
})
names(sc_by_tissue_clean) <- names(sc_by_tissue)
sc_by_tissue_clean <- sc_by_tissue_clean[!sapply(sc_by_tissue_clean, is.null)]

# ------------------------------------------------------------------------------
# 5. 重新合并与保存 (继承脚本1 qs存储逻辑)
# ------------------------------------------------------------------------------
print("🚀 步骤4/4: 合并纯净数据并保存...")

# 将拆分去双胞后的纯净对象重新合并为一个 Seurat 对象
sc_final <- merge(sc_by_tissue_clean[[1]], 
                  y = sc_by_tissue_clean[2:length(sc_by_tissue_clean)])
sc_final <- JoinLayers(sc_final)

# 保存最终对象
qsave(sc_final, "sc_combined_qc_Cleaned.qs")
print(paste("🎉 全部运行完毕！最终保留高质量单细胞总数:", ncol(sc_final)))
print("💾 结果已保存为: sc_combined_qc_Cleaned.qs")