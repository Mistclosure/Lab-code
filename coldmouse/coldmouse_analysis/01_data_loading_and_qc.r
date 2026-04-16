# ==============================================================================
# Script 1: 数据加载与严格质控
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(Matrix)   
library(qs)

# ------------------------------------------------------------------------------
# 1. 设置路径与元数据映射
# ------------------------------------------------------------------------------
data_dir <- "/mnt/disk1/qiuzerui/expriments/coldmouse/sccoldmouse" 
sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")
tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")
group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")

# ------------------------------------------------------------------------------
# 2. 数据读取与 Seurat 对象创建
# ------------------------------------------------------------------------------
sc_list <- list()
full_folder_names <- list.dirs(data_dir, full.names = FALSE, recursive = FALSE)

print("🚀 步骤1/10: 开始读取 10X 矩阵数据并创建对象...")

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
    sc_obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3)
    
    if (ncol(sc_obj) == 0) next 
    
    prefix <- substr(sample, 1, 1); suffix <- substr(sample, 2, 2)
    sc_obj$Tissue <- as.character(tissue_map[prefix])
    sc_obj$Group  <- as.character(group_map[suffix])
    sc_obj$Condition <- ifelse(suffix == "2", "Cold", "NonCold")
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-")
    
    sc_list[[sample]] <- sc_obj
    print(paste("   ✅ 样本入库完成:", sample))
    
  }, error = function(e) { 
    message(paste("   ❌ 处理样本", sample, "时出错:", e$message)) 
  })
}

sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = sample_ids)
sc_combined <- JoinLayers(sc_combined)

# ------------------------------------------------------------------------------
# 3. 严格质控 (QC)
# ------------------------------------------------------------------------------
print("🚀 步骤2/10: 正在进行严格质控...")

sc_combined <- subset(sc_combined, subset = nFeature_RNA >= 100 & nFeature_RNA <= 2500 & percent.mt < 10)

counts_mat <- GetAssayData(sc_combined, layer = "counts")
cells_keep <- colnames(sc_combined)

if ("Pf4" %in% rownames(counts_mat)) {
  cells_keep <- cells_keep[counts_mat["Pf4", cells_keep] == 0]
}
for (rbc_gene in c("Hbb-bs", "Hba-a1")) {
  if (rbc_gene %in% rownames(counts_mat)) {
    cells_keep <- cells_keep[counts_mat[rbc_gene, cells_keep] <= 1]
  }
}
sc_combined <- subset(sc_combined, cells = cells_keep)

ribo_genes <- grep("^Rp[sl]", rownames(sc_combined), value = TRUE, ignore.case = TRUE)
non_ribo_genes <- rownames(sc_combined)[!rownames(sc_combined) %in% ribo_genes]
sc_combined <- subset(sc_combined, features = non_ribo_genes)

print(paste("✅ 严格过滤后剩余细胞数:", ncol(sc_combined)))

# 【新增】保存质控后的对象供下一脚本使用
qsave(sc_combined, "sc_combined_qc.qs")
print("✅ 脚本1 运行完毕，结果已保存为 sc_combined_qc.qs")