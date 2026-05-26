# ==============================================================================
# 适配 Kat8 P60 WT vs KO 数据的完整流程
# 特性：自动处理生物学重复整合 (Harmony) + 智能文件名解析
# ==============================================================================
set.seed(42)

# ------------------------------------------------------------------------------
# 0. 加载必要的 R 包
# ------------------------------------------------------------------------------
library(Seurat)
library(tidyverse)
library(patchwork)
library(openxlsx) 
library(Matrix)    
library(scales)
library(HGNChelper) 
library(dplyr)
library(harmony) # 增加：加载 harmony 包
library(qs)      # 确保加载 qs 以支持后面的 qsave

# --- 加载去双胞专用包 ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("scDblFinder", quietly = TRUE)) BiocManager::install("scDblFinder")
if (!require("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment")

library(scDblFinder)
library(SingleCellExperiment)

# ------------------------------------------------------------------------------
# 1. 设置路径与样本信息
# ------------------------------------------------------------------------------
data_dir <- "/mnt/disk1/qiuzerui/expriments/cortex_Y90C/filtered_counts" 
setwd(data_dir)
# 根据文件夹名定义
sample_folders <- c("kat8-P60-WT-1", "kat8-P60-WT-2", "kat8-P60-Y90C-KO-1", "kat8-P60-Y90C-KO-2")

# ------------------------------------------------------------------------------
# 加载 ScType 核心函数 (本地化版本)
# ------------------------------------------------------------------------------
# 定义本地文件名
sctype_files <- c("gene_sets_prepare.R", "sctype_score_.R")

for (f in sctype_files) {
  local_path <- file.path(data_dir, f)
  
  # 如果本地不存在该文件，则尝试下载一次
  if (!file.exists(local_path)) {
    message(paste("🌐 正在从 GitHub 下载:", f, "..."))
    url <- paste0("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/", f)
    tryCatch({
      download.file(url, destfile = local_path, mode = "wb")
    }, error = function(e) {
      stop(paste("❌ 下载失败，请手动下载文件并放入目录:", data_dir))
    })
  }
  
  # 从本地加载
  source(local_path)
}
print("✅ ScType 核心函数已从本地加载完毕")

# ------------------------------------------------------------------------------
# 2. 读取数据、解析文件名元数据
# ------------------------------------------------------------------------------
sc_list <- list()

print("🚀 步骤1/6: 开始读取数据并解析实验分组...")

for (folder in sample_folders) {
  full_path <- file.path(data_dir, folder)
  
  if(!dir.exists(full_path)) {
    message(paste("⚠️ 跳过：未找到文件夹", folder)); next
  }
  
  print(paste("    正在处理:", folder))
  
  tryCatch({
    counts <- Read10X(data.dir = full_path, gene.column = 1)
    
    if (is.list(counts) && !is(counts, "dgCMatrix")) {
      counts <- if ("Gene Expression" %in% names(counts)) counts$`Gene Expression` else counts[[1]]
    }
    
    sc_obj <- CreateSeuratObject(counts = counts, project = folder, min.cells = 3, min.features = 200)
    
    # 智能解析分组
    if (grepl("WT", folder)) {
      group <- "WT"
    } else if (grepl("KO", folder)) {
      group <- "KO" 
    } else {
      group <- "Unknown"
    }
    
    rep_id <- substr(folder, nchar(folder), nchar(folder))
    
    sc_obj$Orig_Folder <- folder       
    sc_obj$Group       <- group        
    sc_obj$Replicate   <- rep_id       
    sc_obj$SampleID    <- paste0(group, "_Rep", rep_id) 
    
    sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^mt-") 
    
    sc_list[[folder]] <- sc_obj
    print(paste("    ✅ 成功入库: Group=", group, "| Rep=", rep_id, "| Cells:", ncol(sc_obj)))
    
  }, error = function(e) {
    message(paste("    ❌ 出错:", folder, e$message))
  })
}

# 合并所有样本
if (length(sc_list) > 0) {
  print("正在合并样本...")
  sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = names(sc_list))
} else {
  stop("❌ 未读取到数据")
}
sc_combined <- JoinLayers(sc_combined)

# ------------------------------------------------------------------------------
# 3. 统一质控 (QC)
# ------------------------------------------------------------------------------
print("🚀 步骤2/6: 正在进行质控过滤...")
sc_combined <- subset(sc_combined, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 15)

# ------------------------------------------------------------------------------
# 4. 去除双细胞
# ------------------------------------------------------------------------------
print("🚀 步骤3/6: 识别并去除双细胞 (scDblFinder)...")

sce <- as.SingleCellExperiment(sc_combined)
sce <- scDblFinder(sce, samples = "Orig_Folder") 

sc_combined$scDblFinder_class <- sce$scDblFinder.class
n_dbl <- sum(sc_combined$scDblFinder_class == "doublet")
print(paste("    [Result] 总计发现双细胞:", n_dbl, "个"))

sc_combined <- subset(sc_combined, subset = scDblFinder_class == "singlet")

# ------------------------------------------------------------------------------
# 5. 标准流程：降维与聚类 (使用 Harmony 整合)
# ------------------------------------------------------------------------------
print("🚀 步骤4/6: 标准化与生物学重复整合 (Harmony Workflow)...")

obj <- sc_combined
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000)
obj <- ScaleData(obj)
obj <- RunPCA(obj, verbose = FALSE)

# --- 增加 Harmony 整合 ---
print("    正在运行 Harmony 去除批次效应...")
obj <- RunHarmony(obj, group.by.vars = "Orig_Folder")
reduction_to_use <- "harmony" # 后续使用 harmony 的降维结果

obj <- RunUMAP(obj, reduction = reduction_to_use, dims = 1:20)
obj <- FindNeighbors(obj, reduction = reduction_to_use, dims = 1:20)
obj <- FindClusters(obj, resolution = 0.5)

# ------------------------------------------------------------------------------
# 6. ScType 自动注释
# ------------------------------------------------------------------------------
print("🚀 步骤5/6: 运行 ScType 细胞注释...")
db_file_path <- file.path(data_dir, "ScTypeDB_full.xlsx") 

if (file.exists(db_file_path)) {
  gs_list_immune <- gene_sets_prepare(db_file_path, "Brain") 
  
  es.max <- sctype_score(scRNAseqData = as.matrix(GetAssayData(obj, layer="scale.data")), scaled = TRUE, 
                         gs = gs_list_immune$gs_positive, gs2 = gs_list_immune$gs_negative)
  
  cL_resutls <- do.call("rbind", lapply(unique(obj@meta.data$seurat_clusters), function(cl){
    cells_in_cluster <- rownames(obj@meta.data[obj@meta.data$seurat_clusters == cl, ])
    es.max_subset <- es.max[ , cells_in_cluster, drop = FALSE]
    es.max.cl = sort(rowSums(es.max_subset), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl), 10)
  }))
  
  sctype_scores <- cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  
  obj@meta.data$cell_type <- ""
  for(j in unique(sctype_scores$cluster)){
    cl_type <- sctype_scores[sctype_scores$cluster==j, "type"]
    obj@meta.data$cell_type[obj@meta.data$seurat_clusters == j] <- as.character(cl_type)
  }
} else {
  message("⚠️ 未找到 ScType 数据库，跳过注释，使用 Cluster ID 绘图")
  obj$cell_type <- obj$seurat_clusters
}

# 修正特定细胞类型
obj[["cell_type"]][obj[["cell_type"]] == "Tanycytes"] <- "Mature neurons"

# 保存带有 harmony 后缀的 qs 对象
qsave(obj, 'Y90C_harmony.qs')

# ------------------------------------------------------------------------------
# 7. 结果可视化与输出
# ------------------------------------------------------------------------------
print("🚀 步骤6/6: 正在生成 PNG 图片...")

plot_group <- "cell_type"
# 更改出图文件夹名字以区分
plot_dir <- file.path(data_dir, "Results_Plots_Harmony") 
if (!dir.exists(plot_dir)) dir.create(plot_dir)

my_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
  legend.position = "right",
  legend.text = element_text(size = 10),
  legend.title = element_blank()
) 
my_guide <- guides(color = guide_legend(ncol = 1, override.aes = list(size = 4)))

# 1. 总图
p_total <- DimPlot(obj, reduction = "umap", group.by = plot_group, label = TRUE, repel = TRUE) + 
  ggtitle(paste0("Total Merged (Cells: ", ncol(obj), ")")) + my_theme + my_guide

# 2. WT 独立图
obj_wt <- subset(obj, subset = Group == "WT")
p_wt <- DimPlot(obj_wt, reduction = "umap", group.by = plot_group, label = TRUE, repel = TRUE) +
  ggtitle("WT Group") + my_theme + my_guide

# 3. KO 独立图
obj_ko <- subset(obj, subset = Group == "KO")
p_ko <- DimPlot(obj_ko, reduction = "umap", group.by = plot_group, label = TRUE, repel = TRUE) +
  ggtitle("KO Group") + my_theme + my_guide

# 4. 对比图
p_split <- DimPlot(obj, reduction = "umap", group.by = plot_group, split.by = "Group", label = TRUE, ncol = 2) +
  ggtitle("Condition Comparison: WT vs KO") + theme(legend.position = "right") + my_guide

# 所有图片文件名加上 _harmony 后缀
ggsave(file.path(plot_dir, "01_UMAP_Total_harmony.png"), plot = p_total, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "02_UMAP_WT_harmony.png"), plot = p_wt, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "03_UMAP_KO_harmony.png"), plot = p_ko, width = 14, height = 9, dpi = 300, bg = "white")
ggsave(file.path(plot_dir, "04_UMAP_Split_harmony.png"), plot = p_split, width = 20, height = 8, dpi = 300, bg = "white")

# ------------------------------------------------------------------------------
# 8. 细胞比例分析
# ------------------------------------------------------------------------------
print("🚀 步骤7/8: 正在进行细胞比例统计分析...")
stats_dir <- file.path(data_dir, "Results_Stats")
if (!dir.exists(stats_dir)) dir.create(stats_dir)

prop_data <- obj@meta.data %>%
  group_by(Group, cell_type) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  group_by(Group) %>%
  mutate(Total = sum(Count), Percent = Count / Total,
         Label = ifelse(Percent > 0.03, paste0(round(Percent * 100, 1), "%"), ""))

p_barplot <- ggplot(prop_data, aes(x = Group, y = Percent, fill = cell_type)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  geom_text(aes(label = Label), position = position_fill(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Cell Type Proportion: WT vs KO (Harmony)") + theme_classic()

# 表格和图片均加上 _harmony 后缀
ggsave(file.path(stats_dir, "05_Cell_Proportion_Barplot_harmony.png"), plot = p_barplot, width = 8, height = 6)
write.csv(prop_data, file.path(stats_dir, "05_Cell_Proportion_Table_harmony.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 9. 差异分析 (KO vs WT)
# ------------------------------------------------------------------------------
print("🚀 步骤8/8: 正在进行差异基因分析 (KO vs WT)...")
Idents(obj) <- "cell_type"
all_cell_types <- unique(obj$cell_type)
deg_list <- list()

for (ctype in all_cell_types) {
  sub_obj <- subset(obj, idents = ctype)
  Idents(sub_obj) <- "Group"
  group_counts <- table(sub_obj$Group)
  
  if (all(c("WT", "KO") %in% names(group_counts)) && min(group_counts) >= 3) {
    tryCatch({
      markers <- FindMarkers(sub_obj, ident.1 = "KO", ident.2 = "WT", logfc.threshold = 0.25, min.pct = 0.1)
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$cell_type <- ctype
        deg_list[[ctype]] <- markers
      }
    }, error = function(e) { message(paste("误差:", ctype, e$message)) })
  }
}

if (length(deg_list) > 0) {
  all_degs <- do.call(rbind, deg_list)
  wb <- createWorkbook()
  addWorksheet(wb, "All_Combined"); writeData(wb, "All_Combined", all_degs)
  for (ctype_name in names(deg_list)) {
    clean_name <- substr(gsub("[^[:alnum:]]", "_", ctype_name), 1, 30)
    addWorksheet(wb, clean_name); writeData(wb, clean_name, deg_list[[ctype_name]])
  }
  # Excel 输出加上 _harmony 后缀
  saveWorkbook(wb, file.path(stats_dir, "06_DEGs_KO_vs_WT_Full_Report_harmony.xlsx"), overwrite = TRUE)
}

print("🎉 所有分析流程结束！")
k=DimPlot(obj, reduction = "umap", group.by = "Orig_Folder")
ggsave(file.path(stats_dir, "06_source_harmony.png"), plot = k, width = 8, height = 6,dpi = 300)
