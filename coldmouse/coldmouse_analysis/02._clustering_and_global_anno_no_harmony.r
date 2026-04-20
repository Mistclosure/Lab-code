# ==============================================================================
# Script 2: 降维聚类与全局 SingleR 注释 (参数化聚类版) - 承接 Cleaned 数据
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
# 已经去过双胞，注销掉这两个包节约内存
# library(scDblFinder) 
# library(SingleCellExperiment)
library(SingleR)
library(celldex)

# 【关键修改】读取综合质控脚本输出的纯净去双胞数据
sc_combined <- qread("sc_combined_qc_Cleaned.qs")

# ------------------------------------------------------------------------------
# 4. 按组织独立 analysis (Leiden/Louvain)
# ------------------------------------------------------------------------------
print(paste0("🚀 步骤3/10: 使用 ", cluster_method, " 进行降维与聚类 (跳过去双胞，因上游已处理)..."))
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name, method_id) {
  print(paste(">>> 处理组织:", tissue_name, "| 细胞数:", ncol(obj)))
  if(ncol(obj) < 50) return(NULL)
  
  # 【新增】确保 Seurat v5 layer 合并，防止 ScaleData 报错
  obj <- JoinLayers(obj)
  
  # 执行单次完整的 标准化 -> 找高变 -> 缩放 -> PCA
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)
  
  # 聚类
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:50)
  obj <- FindClusters(obj, resolution = 1.0, algorithm = method_id) 
  
  # 【关键修改】删除了冗余的 scDblFinder 和重复的降维代码，直接跑 UMAP
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:50)
  
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

  print(paste("   🔎 正在计算", tissue_name, "的差异表达基因 (基于 Seurat Clusters)..."))
  
  # 【修改点 1】设置 Idents 为聚类 ID 而非 SingleR 标签
  Idents(obj) <- "seurat_clusters"
  
  all_markers <- FindAllMarkers(obj, 
                                only.pos = TRUE, 
                                min.pct = 0.25, 
                                logfc.threshold = 0.25,
                                verbose = FALSE)

  # 【修改点 2】建立 Cluster 与 SingleR 标签的映射关系 (取每个 Cluster 中频数最高的标签)
  cluster_to_label <- as.data.frame(obj@meta.data) %>%
    select(seurat_clusters, SingleR.labels) %>%
    mutate(SingleR.labels = as.character(unlist(SingleR.labels))) %>% # 强制转为标准字符向量，防 list
    group_by(seurat_clusters) %>%
    dplyr::count(SingleR.labels) %>%
    slice_max(n = 1, order_by = n, with_ties = FALSE) %>% # with_ties = FALSE 极其重要，防止平局导致后续 join 行数翻倍
    ungroup() %>%
    select(cluster = seurat_clusters, SingleR_Annotation = SingleR.labels)

  # 【修改点 3】将注释合并到 Marker 结果中
  all_markers <- all_markers %>%
    left_join(cluster_to_label, by = "cluster") %>%
    mutate(cluster_with_anno = paste0("Cluster ", cluster, " (", SingleR_Annotation, ")"))
  
  top20_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)
  
  # 保存 Marker 列表
  write.csv(all_markers, file.path("files", paste0("Markers_All_", cluster_method, "_", tissue_name, ".csv")), row.names = FALSE)
  write.csv(top20_markers, file.path("files", paste0("Markers_Top20_", cluster_method, "_", tissue_name, ".csv")), row.names = FALSE)
  
  print(paste("   ✅", tissue_name, "Marker 列表 (含 SingleR 注释) 已导出。"))

  # 以下绘图部分逻辑保持不变
  all_labels <- sort(unique(obj$SingleR.labels))
  obj$SingleR.labels <- factor(obj$SingleR.labels, levels = all_labels)
  
  sc_by_tissue[[tissue_name]] <- obj
  
  # ==========================================
  # 绘图 1: 基于 SingleR 注释 (按 Annotation 分组)
  # ==========================================
  p_total_anno <- DimPlot(obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "-", cluster_method, "(Annotation)")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())
  
  p_cold_anno <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("Cold_4C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_rt_anno <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("RT_25C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_tn_anno <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
    ggtitle("TN_30C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_final_anno <- (p_total_anno | p_cold_anno) / (p_rt_anno | p_tn_anno) + 
    plot_layout(guides = "collect", axes = "collect") & 
    theme(legend.text = element_text(size = 9))
  
  file_name_anno <- file.path("pictures", paste0("UMAP_Grid_Annotation_", cluster_method, "_", tissue_name, ".png"))
  ggsave(filename = file_name_anno, plot = p_final_anno, width = 15, height = 11, dpi = 300)

  # ==========================================
  # 绘图 2: 基于 Seurat Clusters (按原始聚类分组)
  # ==========================================
  p_total_cls <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "-", cluster_method, "(Clusters)")) + 
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())
  
  p_cold_cls <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "seurat_clusters", label = FALSE) + 
    ggtitle("Cold_4C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_rt_cls <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "seurat_clusters", label = FALSE) + 
    ggtitle("RT_25C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_tn_cls <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "seurat_clusters", label = FALSE) + 
    ggtitle("TN_30C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + 
    NoLegend()
  
  p_final_cls <- (p_total_cls | p_cold_cls) / (p_rt_cls | p_tn_cls) + 
    plot_layout(guides = "collect", axes = "collect") & 
    theme(legend.text = element_text(size = 9))
  
  file_name_cls <- file.path("pictures", paste0("UMAP_Grid_Clusters_", cluster_method, "_", tissue_name, ".png"))
  ggsave(filename = file_name_cls, plot = p_final_cls, width = 15, height = 11, dpi = 300)

  print(paste("   ✅", tissue_name, "的 Annotation 和 Clusters 两版 UMAP 均已保存。"))
}