# ==============================================================================
# 完整流程代码：加载 -> 严格质控 -> 批次校正(Harmony) -> Leiden聚类(Res 1.0) -> 
# 去双胞 -> 提取髓系 -> 基于差异基因(FindAllMarkers)人工校验注释 -> 分组绘图对比
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包 (包含 Leiden 聚类所需的 Python 环境)
# ------------------------------------------------------------------------------
# 指定 Python 环境以运行 Leiden 算法
library(reticulate)
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
library(Seurat)
library(tidyverse)
library(patchwork)
library(Matrix)   
library(qs)
library(leidenbase)
# --- 安装并加载缺失的依赖包 ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")
if (!require("harmony", quietly = TRUE)) install.packages("harmony")
if (!require("SingleR", quietly = TRUE)) BiocManager::install("SingleR")
if (!require("celldex", quietly = TRUE)) BiocManager::install("celldex")
if (!require("leiden", quietly = TRUE)) install.packages("leiden") # 确保加载 leiden 包

library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library(SingleR)
library(celldex)

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

print("🚀 步骤1/6: 开始读取 10X 矩阵数据并创建对象...")

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
print("🚀 步骤2/6: 正在进行严格质控...")

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

# ------------------------------------------------------------------------------
# 4. 按组织独立分析 (Harmony + Leiden)
# ------------------------------------------------------------------------------
print("🚀 步骤3/6: 组织拆分 -> 批次校正 -> Leiden 聚类与去双胞...")
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name) {
  print(paste(">>> 处理组织:", tissue_name, "| 细胞数:", ncol(obj)))
  if(ncol(obj) < 50) return(NULL)
  
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "orig.ident", verbose = FALSE)
  
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:50)
  # 使用 Leiden 聚类 (algorithm = 4), 分辨率 1.0
  obj <- FindClusters(obj, resolution = 1.0, algorithm = 4) 
  
  temp_obj <- JoinLayers(obj)
  sce <- as.SingleCellExperiment(temp_obj)
  sce <- scDblFinder(sce, clusters = TRUE)
  obj$scDblFinder_class <- sce$scDblFinder.class
  rm(temp_obj); gc()
  
  obj <- subset(obj, subset = scDblFinder_class == "singlet")
  
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "orig.ident", verbose = FALSE)
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:50)
  
  return(obj)
}

sc_by_tissue <- lapply(names(sc_by_tissue), function(x) run_standard_pipeline(sc_by_tissue[[x]], x))
names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow")
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]
qsave(sc_by_tissue, 'sc_by_tissue.qs')

# ------------------------------------------------------------------------------
# 5. SingleR 提取并聚类 (加固版)
# ------------------------------------------------------------------------------
print("🚀 步骤4/6: 运行 SingleR 标注并提取髓系细胞...")
sc_by_tissue = qread('sc_by_tissue.qs')
mouse_ref <- celldex::MouseRNAseqData()
myeloid_keywords <- c("Macrophages", "Monocytes", "Granulocytes", "Dendritic cells", "Microglia", "Neutrophils")
sc_myeloid_list <- list()

for (i in seq_along(sc_by_tissue)) {
  tissue_name <- names(sc_by_tissue)[i]
  print(paste("🚀 正在循环处理:", tissue_name))
  obj <- sc_by_tissue[[i]]
  obj <- JoinLayers(obj)
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  pred_res <- SingleR(test = expr_mat, ref = mouse_ref, labels = mouse_ref$label.main)
  obj$SingleR.labels <- pred_res$labels
  sc_by_tissue[[i]] <- obj
  
  is_myeloid <- obj$SingleR.labels %in% myeloid_keywords
  is_myeloid[is.na(is_myeloid)] <- FALSE
  
  if (sum(is_myeloid) >= 20) {
    myeloid_obj <- subset(obj, cells = colnames(obj)[is_myeloid])
    # 髓系子集独立重聚类 (Leiden, Res 1.0)
    myeloid_obj <- FindVariableFeatures(myeloid_obj)
    myeloid_obj <- ScaleData(myeloid_obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
    myeloid_obj <- RunPCA(myeloid_obj, npcs = 50, verbose = FALSE)
    myeloid_obj <- RunHarmony(myeloid_obj, group.by.vars = "orig.ident", verbose = FALSE)
    myeloid_obj <- FindNeighbors(myeloid_obj, reduction = "harmony", dims = 1:30)
    myeloid_obj <- FindClusters(myeloid_obj, resolution = 1.0, algorithm = 4) 
    myeloid_obj <- RunUMAP(myeloid_obj, reduction = "harmony", dims = 1:30)
    sc_myeloid_list[[tissue_name]] <- myeloid_obj
    print(paste("   ✅", tissue_name, "提取到", sum(is_myeloid), "个髓系细胞"))
  }
}
# ------------------------------------------------------------------------------
# 5.5 新增：独立循环导出每个组织每个 Cluster 的 Top 20 Markers
# ------------------------------------------------------------------------------
print("🚀 正在独立导出各组织 Marker 基因表...")

for (tissue in names(sc_myeloid_list)) {
  cat(paste0(">>> 正在提取 [", tissue, "] 的差异基因...\n"))
  
  # 获取对象并确保数据层合并
  obj_tmp <- sc_myeloid_list[[tissue]]
  obj_tmp <- JoinLayers(obj_tmp)
  
  # 计算所有聚类的 Markers
  all_markers <- FindAllMarkers(
    obj_tmp, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  
  # 提取每个 Cluster 的 Top 20
  top20_table <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC) %>%
    as.data.frame()
  
  # 导出 CSV
  output_name <- paste0("Myeloid_Markers_Top20_", tissue, ".csv")
  write.csv(top20_table, output_name, row.names = FALSE)
  
  # 将全量 markers 存入临时列表或直接释放（如果你后面注释还要用 markers，建议在此处保留）
  # 这里为了后面步骤 6 运行顺畅，我们不破坏原有的循环逻辑
  cat(paste0("   ✅ 已保存: ", output_name, "\n"))
}

# ------------------------------------------------------------------------------
# 6. 差异表达分析 (Markers) 与人工校对注释逻辑
# ------------------------------------------------------------------------------
print("🚀 步骤5/6: 运行 FindAllMarkers (Scanpy tl.rank_genes_groups 逻辑) 并注释...")

# 优化后的参考列表：去掉了共有背景基因
refined_myeloid_markers <- list(
  "C1q+ TAM (Resident-like)" = c("C1qa", "C1qb", "C1qc", "Apoe"), 
  "SPP1+ TAM"      = c("Spp1", "Fn1", "Fabp5"),
  "SELENOP+ TAM"   = c("Selenop", "Mrc1", "Cd163"), 
  "EEF1A1+ TAM"    = c("Eef1a1", "Eif3g", "Rps12"), 
  "CX3CR1+ TAM"    = c("Cx3cr1", "Csf1r"),
  "IL1B+ TAM (Inflammatory)" = c("Il1b", "Thbs1", "Ccl3", "Ccl4"),
  "Prolif-TAM"     = c("Mki67", "Top2a", "Stmn1", "Ccnb2"), 
  "IFN-TAM"        = c("Isg15", "Ifit1", "Ifit2", "Cxcl10"), 
  "cDC1"           = c("Clec9a", "Xcr1", "Wdfy4", "Irf8"),
  "cDC2"           = c("Cd209a", "Sirpa", "Irf4", "Klrf1"),
  "mregDC (LAMP3+)"= c("Lamp3", "Ccr7", "Fscn1", "Cd80"), 
  "pDC"            = c("Siglech", "Bst2", "Tlr7", "Irf7"), 
  "Monocytes"      = c("Ly6c2", "Vcan", "F10", "Plac8"),
  "PMN-MDSC/Neutro"= c("S100a8", "S100a9", "Ly6g", "Retnlg"), 
  "M-MDSC"         = c("Arg1", "Ccr2", "Ms4a6c")
)

for (tissue in names(sc_myeloid_list)) {
  obj <- sc_myeloid_list[[tissue]]
  obj <- JoinLayers(obj)
  
  # 寻找差异标记基因
  print(paste(">>> 正在分析组织:", tissue))
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  # 自动化标注逻辑
  annotation_map <- sapply(levels(obj$seurat_clusters), function(cl) {
    top_markers <- markers %>% 
      filter(cluster == cl) %>% 
      slice_max(n = 20, order_by = avg_log2FC) %>% 
      pull(gene)
    
    # 【修复】使用正确的变量名 refined_myeloid_markers
    overlaps <- sapply(refined_myeloid_markers, function(m) length(intersect(top_markers, m)))
    
    if(max(overlaps) > 0) {
      return(names(refined_myeloid_markers)[which.max(overlaps)])
    } else {
      return("Unknown_Myeloid")
    }
  })
  
  # 【修复】使用 unname() 解决 No cell overlap 报错
  obj$Myeloid_Subtype <- unname(annotation_map[as.character(obj$seurat_clusters)])
  
  sc_myeloid_list[[tissue]] <- obj
  print(paste("   ✅", tissue, "精细注释完成。"))
}

qsave(sc_myeloid_list, "scRNA_myeloid_annotated_unbiased.qs")

# ------------------------------------------------------------------------------
# 7. 绘制各组织髓系结果 (2x2 布局)
# ------------------------------------------------------------------------------
print("🎨 步骤6/6: 生成对比绘图...")

for (tissue_name in names(sc_myeloid_list)) {
  print(paste("正在绘制组织:", tissue_name, "..."))
  obj <- sc_myeloid_list[[tissue_name]]
  
  all_types <- sort(unique(obj$Myeloid_Subtype))
  obj$Myeloid_Subtype <- factor(obj$Myeloid_Subtype, levels = all_types)
  
  p_total <- DimPlot(obj, reduction = "umap", group.by = "Myeloid_Subtype", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "- All Myeloid")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_blank())
  
  p_cold <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + 
    ggtitle("Cold_4C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  
  p_rt <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + 
    ggtitle("RT_25C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  
  p_tn <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + 
    ggtitle("TN_30C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  
  p_final <- (p_total | p_cold) / (p_rt | p_tn) + 
    plot_layout(guides = "collect", axes = "collect") & 
    theme(legend.position = "right", legend.text = element_text(size = 9))
  
  file_name <- paste0("Myeloid_Unbiased_Grid_", tissue_name, ".png")
  ggsave(filename = file_name, plot = p_final, width = 15, height = 11, dpi = 300)
}

print("🎉 全部 10X 分析流程运行完毕！")