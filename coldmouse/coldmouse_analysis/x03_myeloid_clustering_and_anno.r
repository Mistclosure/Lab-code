# ==============================================================================
# Script 3: 髓系细胞提取、细分聚类与亚群注释
# 独立运行，依赖数据: sc_by_tissue.qs
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(reticulate)
use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)
library(harmony)

# ------------------------------------------------------------------------------
# 0. 读取前置数据
# ------------------------------------------------------------------------------
print("🚀 正在加载 sc_by_tissue.qs ...")
sc_by_tissue <- qread('sc_by_tissue.qs')

# ------------------------------------------------------------------------------
# 1. 基于 SingleR 结果提取髓系并重新聚类
# ------------------------------------------------------------------------------
print("🚀 步骤6/10: 提取髓系细胞并独立重聚类...")
myeloid_keywords <- c("Macrophages", "Monocytes", "Granulocytes", "Dendritic cells", "Microglia", "Neutrophils")
sc_myeloid_list <- list()

for (i in seq_along(sc_by_tissue)) {
  tissue_name <- names(sc_by_tissue)[i]
  obj <- sc_by_tissue[[i]]
  
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
    print(paste("   ✅", tissue_name, "提取到", sum(is_myeloid), "个髓系细胞并重聚类完毕"))
  }
}

# ------------------------------------------------------------------------------
# 2. 独立循环导出每个组织每个 Cluster 的 Top 20 Markers
# ------------------------------------------------------------------------------
print("🚀 步骤7/10: 正在独立导出各组织 Marker 基因表...")

for (tissue in names(sc_myeloid_list)) {
  cat(paste0(">>> 正在提取 [", tissue, "] 的差异基因...\n"))
  
  obj_tmp <- sc_myeloid_list[[tissue]]
  obj_tmp <- JoinLayers(obj_tmp)
  
  all_markers <- FindAllMarkers(obj_tmp, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  
  top20_table <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log2FC) %>%
    as.data.frame()
  
  output_name <- paste0("Myeloid_Markers_Top20_", tissue, ".csv")
  write.csv(top20_table, output_name, row.names = FALSE)
  cat(paste0("   ✅ 已保存: ", output_name, "\n"))
}

# ------------------------------------------------------------------------------
# 3. 差异表达分析与人工校对注释逻辑
# ------------------------------------------------------------------------------
print("🚀 步骤8/10: 运行注释髓系亚群逻辑...")

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
  
  print(paste(">>> 正在分析组织:", tissue))
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  annotation_map <- sapply(levels(obj$seurat_clusters), function(cl) {
    top_markers <- markers %>% filter(cluster == cl) %>% slice_max(n = 20, order_by = avg_log2FC) %>% pull(gene)
    overlaps <- sapply(refined_myeloid_markers, function(m) length(intersect(top_markers, m)))
    if(max(overlaps) > 0) return(names(refined_myeloid_markers)[which.max(overlaps)]) else return("Unknown_Myeloid")
  })
  
  obj$Myeloid_Subtype <- unname(annotation_map[as.character(obj$seurat_clusters)])
  sc_myeloid_list[[tissue]] <- obj
  print(paste("   ✅", tissue, "精细注释完成。"))
}

qsave(sc_myeloid_list, "scRNA_myeloid_annotated_unbiased.qs")

# ------------------------------------------------------------------------------
# 4. 绘制各组织髓系结果 (2x2 布局)
# ------------------------------------------------------------------------------
print("🎨 步骤9/10: 生成髓系对比绘图...")

for (tissue_name in names(sc_myeloid_list)) {
  print(paste("正在绘制组织:", tissue_name, "..."))
  obj <- sc_myeloid_list[[tissue_name]]
  
  all_types <- sort(unique(obj$Myeloid_Subtype))
  obj$Myeloid_Subtype <- factor(obj$Myeloid_Subtype, levels = all_types)
  
  p_total <- DimPlot(obj, reduction = "umap", group.by = "Myeloid_Subtype", label = TRUE, repel = TRUE) + 
    ggtitle(paste(tissue_name, "- All Myeloid")) + theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())
  
  p_cold <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + ggtitle("Cold_4C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  p_rt <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + ggtitle("RT_25C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  p_tn <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "Myeloid_Subtype", label = FALSE) + ggtitle("TN_30C") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank())
  
  p_final <- (p_total | p_cold) / (p_rt | p_tn) + plot_layout(guides = "collect", axes = "collect") & theme(legend.position = "right", legend.text = element_text(size = 9))
  ggsave(filename = paste0("Myeloid_Unbiased_Grid_", tissue_name, ".png"), plot = p_final, width = 15, height = 11, dpi = 300)
}
print("🎉 脚本 3 运行完毕！")