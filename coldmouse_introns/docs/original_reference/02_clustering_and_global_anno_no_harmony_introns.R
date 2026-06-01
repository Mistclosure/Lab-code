# ==============================================================================
# coldmouse_introns：无 Harmony 的降维聚类、SingleR 全局注释和 marker 输出
# ==============================================================================
#
# 输入路径：
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns
#
# 输入文件：
#   coldmouse_introns_sc_combined_qc_Cleaned.qs
#
# 输出路径：
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns/files
#   /mnt/disk1/qiuzerui/expriments/coldmouse_introns/pictures
#
# 参考脚本：
#   /home/zerui/code/coldmouse/coldmouse_analysis/using/general_preprocessing/clustering_and_global_anno_no_harmony.r
#
# 主要修改点：
#   1. 承接 introns 版本质控输出，而不是原 coldmouse 的 sc_combined_qc_Cleaned.qs。
#   2. 工作目录改为 /mnt/disk1/qiuzerui/expriments/coldmouse_introns。
#   3. 输出对象、marker 表和 UMAP 图片文件名加入 introns 标识，避免覆盖原项目结果。
#   4. 保留原 no-harmony 脚本的 Seurat 预处理、Louvain/Leiden 聚类、
#      SingleR 注释、marker 计算和绘图逻辑。
# ==============================================================================

setwd("/mnt/disk1/qiuzerui/expriments/coldmouse_introns")

# 创建输出文件夹。沿用参考脚本的 pictures/files 输出结构。
dir.create("pictures", showWarnings = FALSE)
dir.create("files", showWarnings = FALSE)

# 聚类方法参数：可在 "leiden" 和 "louvain" 之间切换。
cluster_method <- "louvain"
algo_id <- ifelse(cluster_method == "leiden", 4, 1) # 4 为 Leiden，1 为 Louvain

library(reticulate)
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(leidenbase)
# 上游 introns 质控脚本已经去过双胞，这里不再加载 scDblFinder/SingleCellExperiment。
# library(scDblFinder)
# library(SingleCellExperiment)
library(SingleR)
library(celldex)

# 读取 introns 版本综合质控脚本输出的纯净去双胞数据。
qc_file <- "coldmouse_introns_sc_combined_qc_Cleaned.qs"
if (!file.exists(qc_file)) {
  stop("未找到 introns 质控输出文件：", file.path(getwd(), qc_file))
}
sc_combined <- qread(qc_file)

# ------------------------------------------------------------------------------
# 1. 按组织独立进行标准化、降维和聚类，不运行 Harmony
# ------------------------------------------------------------------------------
print(paste0("步骤1/4：使用 ", cluster_method, " 进行无 Harmony 降维与聚类。"))
sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

run_standard_pipeline <- function(obj, tissue_name, method_id) {
  print(paste(">>> 处理组织:", tissue_name, "| 细胞数:", ncol(obj)))
  if (ncol(obj) < 50) return(NULL)

  # 确保 Seurat v5 layer 合并，避免 ScaleData/GetAssayData 报错。
  obj <- JoinLayers(obj)

  # 标准 Seurat 流程：标准化 -> 高变基因 -> scale -> PCA。
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj, vars.to.regress = c("nCount_RNA", "percent.mt"))
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

  # 按参考脚本参数聚类，不进行 Harmony 校正。
  obj <- FindNeighbors(obj, reduction = "pca", dims = 1:50)
  obj <- FindClusters(obj, resolution = 1.0, algorithm = method_id)
  obj <- RunUMAP(obj, reduction = "pca", dims = 1:50)

  return(obj)
}

sc_by_tissue <- lapply(names(sc_by_tissue), function(x) run_standard_pipeline(sc_by_tissue[[x]], x, algo_id))
names(sc_by_tissue) <- c("Aorta", "PBMC", "BoneMarrow")
sc_by_tissue <- sc_by_tissue[!sapply(sc_by_tissue, is.null)]

if (length(sc_by_tissue) == 0) {
  stop("没有可用于后续分析的组织对象。")
}

# ------------------------------------------------------------------------------
# 2. 运行 SingleR 对所有细胞做全局注释
# ------------------------------------------------------------------------------
print("步骤2/4：运行 SingleR 对所有细胞进行全局注释。")
mouse_ref <- celldex::MouseRNAseqData()

for (i in seq_along(sc_by_tissue)) {
  tissue_name <- names(sc_by_tissue)[i]
  print(paste(">>> 正在运行 SingleR 注释:", tissue_name))
  obj <- sc_by_tissue[[i]]
  obj <- JoinLayers(obj)
  expr_mat <- GetAssayData(obj, assay = "RNA", layer = "data")
  pred_res <- SingleR(test = expr_mat, ref = mouse_ref, labels = mouse_ref$label.main)
  obj$SingleR.labels <- pred_res$labels
  sc_by_tissue[[i]] <- obj
  print(paste("   完成:", tissue_name, "全局注释"))
}

# 保存 introns 版本按组织对象。
sc_by_tissue_file <- paste0("coldmouse_introns_sc_by_tissue_no_harmony_", cluster_method, ".qs")
qsave(sc_by_tissue, sc_by_tissue_file)

# ------------------------------------------------------------------------------
# 3. 绘制 UMAP 并导出 marker 表
# ------------------------------------------------------------------------------
print("步骤3/4：绘制 UMAP 并导出 marker 表。")
sc_by_tissue <- qread(sc_by_tissue_file)

for (tissue_name in names(sc_by_tissue)) {
  print(paste(">>> 正在处理组织:", tissue_name))
  obj <- sc_by_tissue[[tissue_name]]

  min_cells_threshold <- 20
  cell_counts <- table(obj$SingleR.labels)
  keep_labels <- names(cell_counts[cell_counts >= min_cells_threshold])

  if (length(keep_labels) == 0) {
    print(paste("   ", tissue_name, "中没有符合条件的亚群，跳过。"))
    next
  }

  obj <- subset(obj, subset = SingleR.labels %in% keep_labels)

  print(paste("   正在计算", tissue_name, "的差异表达基因，分组依据为 Seurat clusters。"))

  # 设置 Idents 为聚类 ID，而不是 SingleR 标签。
  Idents(obj) <- "seurat_clusters"

  all_markers <- FindAllMarkers(obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                verbose = FALSE)

  # 建立 Cluster 与 SingleR 标签的映射关系，取每个 Cluster 中频数最高的标签。
  cluster_to_label <- as.data.frame(obj@meta.data) %>%
    select(seurat_clusters, SingleR.labels) %>%
    mutate(SingleR.labels = as.character(unlist(SingleR.labels))) %>%
    group_by(seurat_clusters) %>%
    dplyr::count(SingleR.labels) %>%
    slice_max(n = 1, order_by = n, with_ties = FALSE) %>%
    ungroup() %>%
    select(cluster = seurat_clusters, SingleR_Annotation = SingleR.labels)

  # 将注释信息合并回 marker 表，方便人工核对。
  all_markers <- all_markers %>%
    left_join(cluster_to_label, by = "cluster") %>%
    mutate(cluster_with_anno = paste0("Cluster ", cluster, " (", SingleR_Annotation, ")"))

  top20_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC)

  # 保存 introns 版本 marker 表。
  write.csv(all_markers,
            file.path("files", paste0("Introns_Markers_All_", cluster_method, "_", tissue_name, ".csv")),
            row.names = FALSE)
  write.csv(top20_markers,
            file.path("files", paste0("Introns_Markers_Top20_", cluster_method, "_", tissue_name, ".csv")),
            row.names = FALSE)

  print(paste("   完成:", tissue_name, "marker 表已导出。"))

  # 固定 SingleR 标签顺序，保证同一组织内图例顺序稳定。
  all_labels <- sort(unique(obj$SingleR.labels))
  obj$SingleR.labels <- factor(obj$SingleR.labels, levels = all_labels)

  sc_by_tissue[[tissue_name]] <- obj

  # 图 1：按 SingleR 注释分组。
  p_total_anno <- DimPlot(obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) +
    ggtitle(paste(tissue_name, "-", cluster_method, "(Introns Annotation)")) +
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

  file_name_anno <- file.path("pictures", paste0("Introns_UMAP_Grid_Annotation_", cluster_method, "_", tissue_name, ".png"))
  ggsave(filename = file_name_anno, plot = p_final_anno, width = 15, height = 11, dpi = 300)

  # 图 2：按 Seurat clusters 分组。
  p_total_cls <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle(paste(tissue_name, "-", cluster_method, "(Introns Clusters)")) +
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

  file_name_cls <- file.path("pictures", paste0("Introns_UMAP_Grid_Clusters_", cluster_method, "_", tissue_name, ".png"))
  ggsave(filename = file_name_cls, plot = p_final_cls, width = 15, height = 11, dpi = 300)

  print(paste("   完成:", tissue_name, "Annotation 和 Clusters 两版 UMAP 已保存。"))
}

# ------------------------------------------------------------------------------
# 4. 保存带筛选后对象的 introns 版本组织列表
# ------------------------------------------------------------------------------
qsave(sc_by_tissue, paste0("coldmouse_introns_sc_by_tissue_no_harmony_", cluster_method, "_with_markers_input.qs"))
print("步骤4/4：脚本结束。请在 files/ 和 pictures/ 中查看 introns 版本输出。")
