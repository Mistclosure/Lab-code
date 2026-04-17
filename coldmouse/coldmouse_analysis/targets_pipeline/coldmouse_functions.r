# ==============================================================================
# R/functions.R
# 包含所有单细胞分析的 Pipeline 函数
# ==============================================================================

# --- Script 1: 数据加载与严格质控 ---
run_data_loading_and_qc <- function(data_dir) {
  sample_ids <- c("A1", "A2", "A3", "B1", "B2", "B3", "M1", "M2", "M3")
  tissue_map <- c("A" = "Aorta", "B" = "PBMC", "M" = "BoneMarrow")
  group_map <- c("1" = "RT_25C", "2" = "Cold_4C", "3" = "TN_30C")
  
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
  print("✅ 数据加载与质控完成")
  
  return(sc_combined)
}

# --- Script 2: 降维聚类、去双胞与全局 SingleR 注释 ---
run_clustering_and_global_anno <- function(sc_combined) {
  use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
  
  print("🚀 步骤3/10: 组织拆分 -> 批次校正 -> Leiden 聚类与去双胞...")
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
  
  print("🚀 步骤5/10: 绘制所有细胞的聚类 UMAP 并导出 Top Markers...")
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
    
    all_markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
    
    top20_markers <- all_markers %>%
      group_by(cluster) %>%
      slice_max(n = 20, order_by = avg_log2FC)
    
    write.csv(all_markers, paste0("Markers_All_", tissue_name, ".csv"), row.names = FALSE)
    write.csv(top20_markers, paste0("Markers_Top20_", tissue_name, ".csv"), row.names = FALSE)
    
    print(paste("   ✅", tissue_name, "Marker 列表已导出。"))
    
    all_labels <- sort(unique(obj$SingleR.labels))
    obj$SingleR.labels <- factor(obj$SingleR.labels, levels = all_labels)
    
    p_total <- DimPlot(obj, reduction = "umap", group.by = "SingleR.labels", label = TRUE, repel = TRUE) + 
      ggtitle(paste(tissue_name, "- All Cells")) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"), axis.title = element_blank())
    
    p_cold <- DimPlot(subset(obj, subset = Group == "Cold_4C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
      ggtitle("Cold_4C") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_rt <- DimPlot(subset(obj, subset = Group == "RT_25C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
      ggtitle("RT_25C") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_tn <- DimPlot(subset(obj, subset = Group == "TN_30C"), reduction = "umap", group.by = "SingleR.labels", label = FALSE) + 
      ggtitle("TN_30C") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_final <- (p_total | p_cold) / (p_rt | p_tn) + 
      plot_layout(guides = "collect", axes = "collect") & theme(legend.text = element_text(size = 9))
    
    file_name <- paste0("AllCells_UMAP_Grid_", tissue_name, ".png")
    ggsave(filename = file_name, plot = p_final, width = 15, height = 11, dpi = 300)
  }
  print("✅ 全细胞 UMAP 绘图与 Marker 导出完毕。")
  return(sc_by_tissue)
}

# --- Script 3: 髓系细胞提取、细分聚类与亚群注释 ---
run_myeloid_clustering_and_anno <- function(sc_by_tissue) {
  use_condaenv("py_env", conda = "/home/zerui/miniconda3/bin/conda", required = TRUE)
  
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
  print("🎉 髓系分析完毕！")
  return(sc_myeloid_list)
}

# --- Script 4: Fkbp5+ 单核细胞靶向分析 ---
run_fkbp5_targeted_analysis <- function(sc_by_tissue) {
  print("🚀 步骤10/10: 进行 Fkbp5+ 单核细胞靶向分析与绘图...")
  
  for (tissue_name in names(sc_by_tissue)) {
    print(paste(">>> 正在分析组织:", tissue_name, "的单核细胞 Fkbp5 表达..."))
    obj <- sc_by_tissue[[tissue_name]]
    
    if (!"Monocytes" %in% unique(obj$SingleR.labels)) {
      print(paste("   ⚠️", tissue_name, "中未检测到 Monocytes，跳过该组织。"))
      next
    }
    
    mono_obj <- subset(obj, subset = SingleR.labels == "Monocytes")
    
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
    
    fkbp5_counts <- GetAssayData(mono_obj, assay = "RNA", layer = "counts")[target_gene, ]
    mono_obj$Fkbp5_Status <- ifelse(fkbp5_counts > 0, "Fkbp5+ Monocyte", "Fkbp5- Monocyte")
    mono_obj$Fkbp5_Status <- factor(mono_obj$Fkbp5_Status, levels = c("Fkbp5- Monocyte", "Fkbp5+ Monocyte"))
    
    p_vln <- VlnPlot(mono_obj, features = target_gene, group.by = "Fkbp5_Status", pt.size = 0.5, cols = c("#56B4E9", "#E69F00")) +
      ggtitle(paste(tissue_name, "- Fkbp5 Expression in Monocytes")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    ggsave(filename = paste0("Monocytes_Fkbp5_VlnPlot_", tissue_name, ".png"), plot = p_vln, width = 6, height = 5, dpi = 300)
    
    umap_coords <- Embeddings(mono_obj, reduction = "umap")
    x_lims <- range(umap_coords[, 1])
    y_lims <- range(umap_coords[, 2])
    
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
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_feat_rt <- FeaturePlot(subset(mono_obj, subset = Group == "RT_25C"), features = target_gene, reduction = "umap") + 
      scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
      xlim(x_lims) + ylim(y_lims) +
      ggtitle("RT_25C") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_feat_tn <- FeaturePlot(subset(mono_obj, subset = Group == "TN_30C"), features = target_gene, reduction = "umap") + 
      scale_color_gradientn(colors = c("lightgrey", "red"), limits = c(0, max_expr)) +
      xlim(x_lims) + ylim(y_lims) +
      ggtitle("TN_30C") + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_blank()) + NoLegend()
    
    p_feat_final <- (p_feat_total | p_feat_cold) / (p_feat_rt | p_feat_tn) + 
      plot_layout(guides = "collect", axes = "collect") 
    
    ggsave(filename = paste0("Monocytes_Fkbp5_FeaturePlot_Grid_", tissue_name, ".png"), plot = p_feat_final, width = 15, height = 11, dpi = 300)
    
    print(paste("   ✅", tissue_name, "单核细胞 Fkbp5 绘图完成！"))
  }
  return("Fkbp5_Analysis_Done")
}

# --- Script 5: PBMC Cold_4C 单核细胞组内异质性分析 ---
run_pbmc_dea_enrichment <- function(sc_by_tissue) {
  if (!"PBMC" %in% names(sc_by_tissue)) stop("❌ 未找到 PBMC 数据！")
  
  pbmc_obj <- sc_by_tissue[["PBMC"]]
  pbmc_obj <- JoinLayers(pbmc_obj)
  
  print("🚀 正在提取 PBMC 中 Cold_4C 组的单核细胞...")
  target_cells <- subset(pbmc_obj, subset = Group == "Cold_4C" & SingleR.labels == "Monocytes")
  
  if (ncol(target_cells) < 10) stop("❌ 符合条件的细胞数太少，无法进行差异分析！")
  
  target_gene <- grep("Fkbp5", rownames(target_cells), ignore.case = TRUE, value = TRUE)[1]
  fkbp5_counts <- GetAssayData(target_cells, assay = "RNA", layer = "counts")[target_gene, ]
  
  target_cells$Fkbp5_Status <- ifelse(fkbp5_counts > 0, "Fkbp5+", "Fkbp5-")
  
  print("✅ 分组完成，细胞数量统计：")
  print(table(target_cells$Fkbp5_Status))
  
  Idents(target_cells) <- "Fkbp5_Status"
  
  print("🚀 正在计算差异表达基因 (Fkbp5+ vs Fkbp5-)...")
  deg_results <- FindMarkers(
    target_cells, 
    ident.1 = "Fkbp5+", 
    ident.2 = "Fkbp5-", 
    logfc.threshold = 0.25,
    min.pct = 0.05
  )
  
  deg_results <- deg_results %>% 
    rownames_to_column(var = "gene") %>%
    arrange(desc(avg_log2FC))
  
  write.csv(deg_results, "PBMC_Cold4C_Monocytes_Fkbp5Plus_vs_Minus_DEGs.csv", row.names = FALSE)
  print("✅ 差异基因计算完成，已导出 CSV！")
  
  deg_results$Significance <- "Not Sig"
  deg_results$Significance[deg_results$avg_log2FC > 0.25 & deg_results$p_val_adj < 0.05] <- "Up in Fkbp5+"
  deg_results$Significance[deg_results$avg_log2FC < -0.25 & deg_results$p_val_adj < 0.05] <- "Up in Fkbp5-"
  
  p_volcano <- ggplot(deg_results, aes(x = avg_log2FC, y = -log10(p_val_adj), color = Significance)) +
    geom_point(alpha = 0.8, size = 1.5) +
    scale_color_manual(values = c("Not Sig" = "grey", "Up in Fkbp5+" = "#E64B35", "Up in Fkbp5-" = "#4DBBD5")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    theme_classic() +
    ggtitle("Fkbp5+ vs Fkbp5- Monocytes in Cold_4C PBMC") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave("PBMC_Cold4C_Monocytes_Fkbp5_Volcano.png", plot = p_volcano, width = 7, height = 6, dpi = 300)
  
  print("🚀 正在进行 GO (BP) 富集分析...")
  genes_up <- deg_results %>% filter(avg_log2FC > 0.25 & p_val_adj < 0.05) %>% pull(gene)
  
  if (length(genes_up) > 10) {
    go_bp_up <- enrichGO(
      gene          = genes_up,
      OrgDb         = org.Mm.eg.db,
      keyType       = 'SYMBOL',
      ont           = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )
    
    if (!is.null(go_bp_up) && nrow(go_bp_up@result) > 0) {
      p_go_up <- dotplot(go_bp_up, showCategory = 15, title = "GO BP Enriched in Fkbp5+ Monocytes")
      ggsave("PBMC_Cold4C_Monocytes_Fkbp5Plus_GO_BP_Up.png", plot = p_go_up, width = 8, height = 7, dpi = 300)
      write.csv(as.data.frame(go_bp_up), "PBMC_Cold4C_Monocytes_Fkbp5Plus_GO_BP_Up.csv", row.names = FALSE)
      print("   ✅ Fkbp5+ 上调基因 GO 富集绘图完成！")
    } else {
      print("   ⚠️ Fkbp5+ 上调基因未富集到显著的 GO 通路。")
    }
  } else {
    print("   ⚠️ Fkbp5+ 上调基因太少 (<10)，跳过富集分析。")
  }
  
  return("DEA_Enrichment_Done")
}