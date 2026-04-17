# ==============================================================================
# Script 2: 降维聚类、去双胞与全局 SingleR 注释 (参数化聚类版)
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 【新增】创建输出文件夹
dir.create("pictures", showWarnings = FALSE)
dir.create("files", showWarnings = FALSE)

# 【参数配置】在此处切换 "leiden" 或 "louvain"
cluster_method <- "louvain" 
algo_id <- ifelse(cluster_method == "leiden", 4, 1) # 4 为 Leiden, 1 为 Louvain

library(reticulate)
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(leidenbase)
library(scDblFinder)
library(SingleCellExperiment)
library(harmony)
library(SingleR)
library(celldex)

# 读取上一脚本的结果
sc_combined <- qread("sc_combined_qc.qs")

# ------------------------------------------------------------------------------
# 4. 按组织独立 analysis (Harmony + Leiden/Louvain)
# ------------------------------------------------------------------------------
print(paste0("🚀 步骤3/10: 使用 ", cluster_method, " 进行聚类与去双胞..."))
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name, method_id) {
  print(paste(">>> 处理组织:", tissue_name, "| 细胞数:", ncol(obj)))
  if(ncol(obj) < 50) return(NULL)
  
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  obj <- RunHarmony(obj, group.by.vars = "orig.ident", verbose = FALSE)
  
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:50)
  
  # 使用传入的 method_id 参数
  obj <- FindClusters(obj, resolution = 1.0, algorithm = method_id) 
  
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

# 将 algo_id 传入函数
sc_by_tissue <- lapply(names(sc_by_tissue), function(x) run_standard_pipeline(sc_by_tissue[[x]], x, algo_id))
names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow")
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]

# ------------------------------------------------------------------------------
# 5. [前置调整] 运行 SingleR 对所有细胞进行初步注释
# ------------------------------------------------------------------------------
print("🚀 步骤4/10: 运行 SingleR 对所有细胞进行全局标注...")
mouse_ref <- celldex::MouseRNAseqData()

for (i in seq_along(sc_by_tissue)) {
  tissue_name <- names(sc_by_tissue)[i]
  print(paste("🚀 正在运行 SingleR 注释:", tissue_name))
  obj <- sc_by_tissue[[i]]
  obj <- JoinLayers(obj)
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  pred_res <- SingleR(test = expr_mat, ref = mouse_ref, labels = mouse_ref$label.main)
  obj$SingleR.labels <- pred_res$labels
  sc_by_tissue[[i]] <- obj
  print(paste("   ✅", tissue_name, "全局注释完成"))
}
# 保存文件名增加聚类方法后缀 (保存在当前目录)
qsave(sc_by_tissue, paste0('sc_by_tissue_', cluster_method, '.qs'))

# ------------------------------------------------------------------------------
# 6. 绘制所有细胞的聚类图片 (分组织、分温度) 与 输出 Top Markers
# ------------------------------------------------------------------------------
print("🚀 步骤5/10: 绘制所有细胞的聚类 UMAP 并导出 Top Markers...")
sc_by_tissue = qread(paste0('sc_by_tissue_', cluster_method, '.qs'))

for (tissue_name in names(sc_by_tissue)) {
  print(paste(">>> 正在处理组织:", tissue_name))
  obj <- sc_by_tissue[[tissue_name]]
  
  min_cells_threshold <- 20
  cell_counts <- table(obj$SingleR.labels)
  keep_labels <- names(cell_counts[cell_counts >= min_cells_threshold])
  
  if (length(keep_labels) == 0) {
    print(paste("   ⚠️", tissue_name, "中没有符合条件的亚群，跳过。"))
    next
  }
  
  obj <- subset(obj, subset = SingleR.labels %in% keep_labels)

  print(paste("   🔎 正在计算", tissue_name, "的差异表达基因..."))
  
  Idents(obj) <- "SingleR.labels"
  
  all_markers <- FindAllMarkers(obj, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25,
                                verbose = FALSE)
  
  top20_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  # 【修改】保存 Marker 列表到 files 文件夹
  write.csv(all_markers, file.path("files", paste0("Markers_All_", cluster_method, "_", tissue_name, ".csv")), row.names = FALSE)
  write.csv(top20_markers, file.path("files", paste0("Markers_Top20_", cluster_method, "_", tissue_name, ".csv")), row.names = FALSE)
  
  print(paste("   ✅", tissue_name, "Marker 列表已导出。"))

  all_labels <- sort(unique(obj$SingleR.labels))
  obj$SingleR.labels <- factor(obj$SingleR.labels, levels = all_labels)
  
  # Total 图
  p_total <- DimPlot(obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "-", cluster_method, "(All)")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title = element_blank())
  
  p_cold <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("Cold_4C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_rt <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("RT_25C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_tn <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("TN_30C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_final <- (p_total | p_cold) / (p_rt | p_tn) + 
    plot_layout(guides = "collect", axes = "collect") & 
    theme(legend.text = element_text(size = 9))
  
  # 【修改】图片输出到 pictures 文件夹
  file_name <- file.path("pictures", paste0("UMAP_Grid_", cluster_method, "_", tissue_name, ".png"))
  ggsave(filename = file_name, plot = p_final, width = 15, height = 11, dpi = 300)
}
print("✅ 全细胞 UMAP 绘图与 Marker 导出完毕。")