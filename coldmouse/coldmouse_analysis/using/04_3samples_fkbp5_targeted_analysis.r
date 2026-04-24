# ==============================================================================
# Script 4: Fkbp5+ 单核细胞靶向分析
# 独立运行，依赖数据: sc_by_tissue.qs
# ==============================================================================
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

library(Seurat)
library(tidyverse)
library(patchwork)
library(qs)

# ------------------------------------------------------------------------------
# 0. 读取前置数据
# ------------------------------------------------------------------------------
print("🚀 正在加载 sc_by_tissue.qs ...")
sc_by_tissue <- qread('sc_by_tissue.qs')

# ------------------------------------------------------------------------------
# 1. Fkbp5+ / Fkbp5- 单核细胞专属分析 (VlnPlot & UMAP FeaturePlot)
# ------------------------------------------------------------------------------
print("🚀 步骤10/10: 进行 Fkbp5+ 单核细胞靶向分析与绘图...")

for (tissue_name in names(sc_by_tissue)) {
  print(paste(">>> 正在分析组织:", tissue_name, "的单核细胞 Fkbp5 表达..."))
  obj <- sc_by_tissue[[tissue_name]]
  
  # 【独立运行关键适配】直接依据全局 SingleR 注释进行 Subset，解耦对脚本文本3的依赖
  if (!"Monocytes" %in% unique(obj$SingleR.labels)) {
    print(paste("   ⚠️", tissue_name, "中未检测到 Monocytes，跳过该组织。"))
    next
  }
  
  # 提取 Monocytes
  mono_obj <- subset(obj, subset = SingleR.labels == "Monocytes")
  
  # 确定 Fkbp5 基因名
  if (!"Fkbp5" %in% rownames(mono_obj)) {
    real_name <- grep("Fkbp5", rownames(mono_obj), ignore.case = TRUE, value = TRUE)
    if(length(real_name) > 0) {
      target_gene <- real_name[1]
    } else {
      print(paste("   ⚠️", tissue_name, "中未检测到 Fkbp5 基因，跳过。"))
      next
    }
  } else {
    target_gene <- "Fkbp5"
  }
  
  # 提取 counts 表达量判定正负
  fkbp5_counts <- GetAssayData(mono_obj, assay = "RNA", layer = "counts")[target_gene, ]
  mono_obj$Fkbp5_Status <- ifelse(fkbp5_counts > 0, "Fkbp5+ Monocyte", "Fkbp5- Monocyte")
  mono_obj$Fkbp5_Status <- factor(mono_obj$Fkbp5_Status, levels = c("Fkbp5- Monocyte", "Fkbp5+ Monocyte"))
  
  # --- (1) 小提琴图 ---
  p_vln <- VlnPlot(mono_obj, features = target_gene, group.by = "Fkbp5_Status", pt.size = 0.5, cols = c("#56B4E9", "#E69F00")) +
    ggtitle(paste(tissue_name, "- Fkbp5 Expression in Monocytes")) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename = paste0("Monocytes_Fkbp5_VlnPlot_", tissue_name, ".png"), plot = p_vln, width = 6, height = 5, dpi = 300)
  
  # ============================================================================
  # --- (2) UMAP FeaturePlot (修复坐标轴与颜色标尺对齐) ---
  # ============================================================================
  
  # 1. 提取全局 UMAP 的坐标极值，固定 X 和 Y 轴 (防止子图放大/变形)
  umap_coords <- Embeddings(mono_obj, reduction = "umap")
  x_lims <- range(umap_coords[, 1])
  y_lims <- range(umap_coords[, 2])
  
  # 2. 提取靶基因的最大表达量，固定颜色标尺 (防止低表达组出现“假红”)
  max_expr <- max(FetchData(mono_obj, vars = target_gene)[, 1])
  
  p_feat_total <- FeaturePlot(mono_obj, features = target_gene, reduction = "umap") + 
    scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
    xlim(x_lims) + ylim(y_lims) +
    ggtitle(paste(tissue_name, "Monocytes - Total")) +
    theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())
  
  p_feat_cold <- FeaturePlot(subset(mono_obj, subset = Group == "Cold_4C"), features = target_gene, reduction = "umap") + 
    scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
    xlim(x_lims) + ylim(y_lims) +
    ggtitle("Cold_4C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()
  
  p_feat_rt <- FeaturePlot(subset(mono_obj, subset = Group == "RT_25C"), features = target_gene, reduction = "umap") + 
    scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
    xlim(x_lims) + ylim(y_lims) +
    ggtitle("RT_25C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()
  
  p_feat_tn <- FeaturePlot(subset(mono_obj, subset = Group == "TN_30C"), features = target_gene, reduction = "umap") + 
    scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
    xlim(x_lims) + ylim(y_lims) +
    ggtitle("TN_30C") + 
    theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) +
    NoLegend()
  
  p_feat_final <- (p_feat_total | p_feat_cold) / (p_feat_rt | p_feat_tn) + 
    plot_layout(guides = "collect", axes = "collect") 
  
  ggsave(filename = paste0("Monocytes_Fkbp5_FeaturePlot_Grid_", tissue_name, ".png"), plot = p_feat_final, width = 15, height = 11, dpi = 300)
  
  print(paste("   ✅", tissue_name, "单核细胞 Fkbp5 绘图完成！"))
}

print("🎉 脚本 4 运行完毕！")