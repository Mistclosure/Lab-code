# ==============================================================================
# hdWGCNA 分析流程：恶性细胞子集 + CiliaHub 基因 + Signed 共表达网络
# 修正版：解决 Assays() / corFnc="pearson" / custom gene_list 等问题
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(Seurat)
  library(tidyverse)
  library(patchwork)
  library(WGCNA)
  library(hdWGCNA)
  library(qs)
})

# ==============================================================================
# 1. 全局参数
# ==============================================================================

WORK_DIR <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"
SEURAT_OBJ_PATH <- file.path(WORK_DIR, "Malignant_RNA_assay.qs")
GENE_LIST_PATH <- "/mnt/disk1/qiuzerui/downloads/signature/ciliahub_genes_list.csv"

setwd(WORK_DIR)

QS_DIR <- "qs"
PLOTS_DIR <- "plots"
FILES_DIR <- "files"

dir.create(QS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FILES_DIR, showWarnings = FALSE, recursive = TRUE)

# hdWGCNA 实验名称
WGCNA_NAME <- "CiliaHub_hdWGCNA_signed"

# Seurat assay 设置
ASSAY <- "RNA"

# 元数据列
# 你的原代码使用 seurat_clusters，因此这里默认使用 seurat_clusters。
CELL_GROUP_COL <- "seurat_clusters"

# 如果你的对象中有样本列，比如 orig.ident / sample / patient，可以填在这里。
# 脚本会自动检测这些列；如果没有，就只按 cluster 构建 metacell。
CANDIDATE_SAMPLE_COLS <- c("orig.ident", "sample", "Sample", "patient", "Patient", "dataset")

# 构建 metacell 的参数
METACELL_K <- 25
METACELL_MAX_SHARED <- 10
METACELL_MIN_CELLS <- 50
METACELL_TARGET <- 1000
METACELL_REDUCTION <- "pca"
METACELL_DIMS <- 1:30

# WGCNA 参数
NETWORK_TYPE <- "signed"
TOM_TYPE <- "signed"

# 注意：
# hdWGCNA / WGCNA 的 corFnc 需要传函数名："cor" 或 "bicor"，
# 不能传 "pearson"，否则会报：没有状态为 "function" 的 "pearson" 目标对象。
# COR_METHOD <- "pearson"
# COR_FNC <- "cor"
# COR_OPTIONS <- "use = 'pairwise.complete.obs'"

# 如果想使用 bicor，可以改成：
# COR_METHOD <- "bicor"
# COR_FNC <- "bicor"
# COR_OPTIONS <- "use='p'"

powers <- c(1:10, seq(12, 30, by = 2))

# 模块提取参数
TARGET_COLORS <- c("brown", "blue", "turquoise", "red", "green")
TARGET_MODULES_DIR <- file.path(FILES_DIR, "Target_Modules_hdWGCNA_单独列表")
dir.create(TARGET_MODULES_DIR, showWarnings = FALSE, recursive = TRUE)

set.seed(9)

# ==============================================================================
# 2. 载入数据
# ==============================================================================

cat("正在加载恶性细胞 Seurat 对象...\n")
malig_seurat <- qread(SEURAT_OBJ_PATH)

# Seurat v5 多 layer 对象需要 JoinLayers；如果不是 v5 或不需要，会自动跳过。
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  malig_seurat <- tryCatch(
    JoinLayers(malig_seurat),
    error = function(e) {
      message("JoinLayers 跳过：", conditionMessage(e))
      malig_seurat
    }
  )
}

# 稳定检查 assay：不要用 ASSAY %in% Assays(obj)，某些环境会触发 match 非向量报错。
available_assays <- names(malig_seurat@assays)
available_assays <- unname(as.character(available_assays))
cat("当前 Seurat 对象包含的 assay：", paste(available_assays, collapse = ", "), "\n")

if (!(ASSAY %in% available_assays)) {
  stop(
    "Seurat 对象中没有 assay: ", ASSAY,
    "\n当前可用 assay: ", paste(available_assays, collapse = ", ")
  )
}
DefaultAssay(malig_seurat) <- ASSAY

if (!(CELL_GROUP_COL %in% colnames(malig_seurat@meta.data))) {
  stop(
    "metadata 中没有列：", CELL_GROUP_COL,
    "\n当前 metadata 列名包括：\n",
    paste(colnames(malig_seurat@meta.data), collapse = ", ")
  )
}

malig_seurat@meta.data[[CELL_GROUP_COL]] <- unname(
  as.character(malig_seurat@meta.data[[CELL_GROUP_COL]])
)

# ==============================================================================
# 3. 自动检查 PCA / Normalize
# ==============================================================================

# hdWGCNA 的 MetacellsByGroups 默认基于 PCA 邻居构建 metacell。
# 如果对象没有 PCA，这里自动补上。
if (!(METACELL_REDUCTION %in% names(malig_seurat@reductions))) {
  cat("未检测到 PCA，正在运行 NormalizeData / FindVariableFeatures / ScaleData / RunPCA...\n")
  malig_seurat <- NormalizeData(malig_seurat, assay = ASSAY, verbose = FALSE)
  malig_seurat <- FindVariableFeatures(malig_seurat, assay = ASSAY, nfeatures = 3000, verbose = FALSE)
  malig_seurat <- ScaleData(malig_seurat, assay = ASSAY, verbose = FALSE)
  malig_seurat <- RunPCA(
    malig_seurat,
    assay = ASSAY,
    npcs = max(METACELL_DIMS),
    verbose = FALSE
  )
} else {
  cat("检测到已有 PCA，跳过 PCA 计算。\n")
}

# ==============================================================================
# 4. 读取 CiliaHub 目标基因
# ==============================================================================

cat("正在读取 CiliaHub 基因列表...\n")
ciliahub_genes <- read.csv(
  GENE_LIST_PATH,
  header = TRUE,
  check.names = FALSE,
  row.names = 1
)

ciliahub_genes_vec <- rownames(ciliahub_genes)
ciliahub_genes_vec <- unname(as.character(ciliahub_genes_vec))

intersect_genes <- intersect(ciliahub_genes_vec, rownames(malig_seurat))
intersect_genes <- unname(as.character(intersect_genes))

cat("CiliaHub 输入基因数：", length(ciliahub_genes_vec), "\n")
cat("Seurat 对象中可用的 CiliaHub 交集基因数：", length(intersect_genes), "\n")

if (length(intersect_genes) < 20) {
  stop("交集基因数少于 20，无法稳定构建共表达网络。请检查基因 ID 是否一致，例如 symbol / Ensembl。")
}

write.csv(
  data.frame(Gene = intersect_genes),
  file = file.path(FILES_DIR, "CiliaHub_intersect_genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

# ==============================================================================
# 5. 设置 hdWGCNA
# ==============================================================================

cat("正在 SetupForWGCNA...\n")

# 使用官方推荐的 custom gene_list 写法。
malig_seurat <- SetupForWGCNA(
  malig_seurat,
  gene_select = "custom",
  gene_list = intersect_genes,
  wgcna_name = WGCNA_NAME
)

# 新建一个整体分组列：所有恶性细胞共同构建一个网络。
# metacell 构建仍可按 cluster / sample 分组，避免不合理混合。
NETWORK_GROUP_COL <- "hdWGCNA_network_group"
NETWORK_GROUP_NAME <- "all_malignant"

malig_seurat@meta.data[[NETWORK_GROUP_COL]] <- NETWORK_GROUP_NAME
malig_seurat@meta.data[[NETWORK_GROUP_COL]] <- unname(
  as.character(malig_seurat@meta.data[[NETWORK_GROUP_COL]])
)

# 自动寻找样本列
sample_col <- CANDIDATE_SAMPLE_COLS[CANDIDATE_SAMPLE_COLS %in% colnames(malig_seurat@meta.data)]
sample_col <- ifelse(length(sample_col) > 0, sample_col[1], NA)

if (!is.na(sample_col)) {
  cat("检测到样本列：", sample_col, "\n")
  malig_seurat@meta.data[[sample_col]] <- unname(as.character(malig_seurat@meta.data[[sample_col]]))
  GROUP_BY_COLS <- c(NETWORK_GROUP_COL, CELL_GROUP_COL, sample_col)
} else {
  cat("未检测到样本列，仅按整体分组 + cluster 构建 metacell。\n")
  GROUP_BY_COLS <- c(NETWORK_GROUP_COL, CELL_GROUP_COL)
}

cat("MetacellsByGroups 使用的 group.by：", paste(GROUP_BY_COLS, collapse = ", "), "\n")

# ==============================================================================
# 6. 构建 metacells 并标准化
# ==============================================================================

cat("正在构建 metacells...\n")
malig_seurat <- MetacellsByGroups(
  seurat_obj = malig_seurat,
  group.by = GROUP_BY_COLS,
  ident.group = CELL_GROUP_COL,
  assay = ASSAY,
  layer = "counts",
  slot = "counts",
  reduction = METACELL_REDUCTION,
  dims = METACELL_DIMS,
  k = METACELL_K,
  max_shared = METACELL_MAX_SHARED,
  min_cells = METACELL_MIN_CELLS,
  target_metacells = METACELL_TARGET,
  wgcna_name = WGCNA_NAME,
  verbose = TRUE
)

cat("正在标准化 metacells...\n")
malig_seurat <- NormalizeMetacells(
  malig_seurat,
  wgcna_name = WGCNA_NAME
)

# ==============================================================================
# 7. 设置表达矩阵
# ==============================================================================

cat("正在 SetDatExpr...\n")
malig_seurat <- SetDatExpr(
  malig_seurat,
  group_name = NETWORK_GROUP_NAME,
  group.by = NETWORK_GROUP_COL,
  assay = ASSAY,
  layer = "data",
  slot = "data",
  wgcna_name = WGCNA_NAME
)

# ==============================================================================
# 8. 软阈值测试
# ==============================================================================

cat("正在测试 soft powers...\n")
malig_seurat <- TestSoftPowers(
  malig_seurat,
  powers = powers,
  networkType = NETWORK_TYPE,
  wgcna_name = WGCNA_NAME
)

# 保存 soft power 图
cat("正在保存 soft power 图...\n")
png(
  filename = file.path(PLOTS_DIR, "hdWGCNA_Signed_SoftPowers.png"),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)
print(PlotSoftPowers(malig_seurat, wgcna_name = WGCNA_NAME))
dev.off()

# 自动选择 soft power：
# ConstructNetwork 默认会根据 TestSoftPowers 的结果自动选择。
# 如果你想固定软阈值，手动设置 soft_power <- 12 即可。
soft_power <- 12

# ==============================================================================
# 9. 构建 signed 共表达网络
# ==============================================================================

cat("正在构建 hdWGCNA signed network...\n")

# 注意：ConstructNetwork 通过 ... 传给 WGCNA::blockwiseConsensusModules。
# 这里 corFnc 必须是 "cor" 或 "bicor"，不能是 "pearson"。
malig_seurat <- ConstructNetwork(
  malig_seurat,
  soft_power = soft_power,
  networkType = NETWORK_TYPE,
  TOMType = TOM_TYPE,
  minModuleSize = MIN_MODULE_SIZE,
  mergeCutHeight = MERGE_CUT_HEIGHT,
  wgcna_name = WGCNA_NAME
)

# 保存 dendrogram
cat("正在保存模块聚类树...\n")
png(
  filename = file.path(PLOTS_DIR, "hdWGCNA_Signed_Module_Dendrogram.png"),
  width = 10,
  height = 6,
  units = "in",
  res = 300
)
PlotDendrogram(malig_seurat, main = "hdWGCNA module dendrogram", wgcna_name = WGCNA_NAME)
dev.off()
# ==============================================================================
# 9.5 WGCNA Network Heatmap Plot
# ==============================================================================

cat("正在绘制 WGCNA Network Heatmap Plot...\n")

# 提取 hdWGCNA 使用的表达矩阵：metacell/sample x gene
datExpr <- GetDatExpr(
  malig_seurat,
  wgcna_name = WGCNA_NAME
)

# 提取模块信息
modules_for_tom <- GetModules(
  malig_seurat,
  wgcna_name = WGCNA_NAME
)

gene_col <- if ("gene_name" %in% colnames(modules_for_tom)) {
  "gene_name"
} else if ("gene" %in% colnames(modules_for_tom)) {
  "gene"
} else {
  colnames(modules_for_tom)[1]
}

color_col <- if ("color" %in% colnames(modules_for_tom)) {
  "color"
} else if ("module" %in% colnames(modules_for_tom)) {
  "module"
} else {
  colnames(modules_for_tom)[2]
}

moduleColors_for_tom <- modules_for_tom[[color_col]]
names(moduleColors_for_tom) <- modules_for_tom[[gene_col]]

# 对齐 datExpr 中的基因顺序
common_genes_tom <- intersect(colnames(datExpr), names(moduleColors_for_tom))
datExpr_tom <- datExpr[, common_genes_tom, drop = FALSE]
moduleColors_tom <- moduleColors_for_tom[common_genes_tom]

# 如果基因太多，TOM 热图会很慢；这里默认最多画 2000 个基因
MAX_TOM_GENES <- 2000

if (ncol(datExpr_tom) > MAX_TOM_GENES) {
  cat("TOM 基因数较多，按模块非 grey 优先抽取前 ", MAX_TOM_GENES, " 个基因用于热图...\n")
  
  non_grey_genes <- names(moduleColors_tom)[moduleColors_tom != "grey"]
  selected_genes <- head(non_grey_genes, MAX_TOM_GENES)
  
  datExpr_tom <- datExpr_tom[, selected_genes, drop = FALSE]
  moduleColors_tom <- moduleColors_tom[selected_genes]
}

# 计算 TOM
adjacency_tom <- adjacency(
  datExpr_tom,
  power = soft_power,
  type = NETWORK_TYPE
)

TOM <- TOMsimilarity(
  adjacency_tom,
  TOMType = TOM_TYPE
)

dissTOM <- 1 - TOM
geneTree_tom <- hclust(as.dist(dissTOM), method = "average")

plotTOM <- dissTOM^4
diag(plotTOM) <- NA

png(
  filename = file.path(PLOTS_DIR, "hdWGCNA_Signed_Network_Heatmap.png"),
  width = 9,
  height = 9,
  units = "in",
  res = 300
)

TOMplot(
  plotTOM,
  geneTree_tom,
  moduleColors_tom,
  main = "Network heatmap plot, all selected genes"
)

dev.off()

qsave(
  list(
    TOM = TOM,
    dissTOM = dissTOM,
    geneTree = geneTree_tom,
    moduleColors = moduleColors_tom
  ),
  file = file.path(QS_DIR, "hdWGCNA_TOM_heatmap_data.qs")
)

cat("WGCNA Network Heatmap Plot 已保存至：",
    file.path(PLOTS_DIR, "hdWGCNA_Signed_Network_Heatmap.png"), "\n")

# ==============================================================================
# 10. 模块特征基因与 kME
# ==============================================================================

cat("正在计算 Module Eigengenes...\n")
malig_seurat <- ModuleEigengenes(
  malig_seurat,
  assay = ASSAY,
  wgcna_name = WGCNA_NAME
)

cat("正在计算 Module Connectivity / kME...\n")
malig_seurat <- ModuleConnectivity(
  malig_seurat,
  group.by = NETWORK_GROUP_COL,
  group_name = NETWORK_GROUP_NAME,
  corFnc = COR_FNC,
  corOptions = COR_OPTIONS,
  assay = ASSAY,
  layer = "data",
  slot = "data",
  wgcna_name = WGCNA_NAME
)

# ==============================================================================
# 11. 结果提取与导出
# ==============================================================================

cat("正在提取模块基因信息...\n")
modules <- GetModules(malig_seurat, wgcna_name = WGCNA_NAME)

# hdWGCNA 的 modules 通常包含 gene_name / module / color / kME 等列。
write.csv(
  modules,
  file = file.path(FILES_DIR, "Module_Genes_hdWGCNA_signed_full.csv"),
  row.names = FALSE,
  quote = FALSE
)

# 兼容你原来下游分析格式：Gene / Module
gene_col <- if ("gene_name" %in% colnames(modules)) "gene_name" else colnames(modules)[1]
color_col <- if ("color" %in% colnames(modules)) "color" else if ("module" %in% colnames(modules)) "module" else colnames(modules)[2]

gene_module_df <- data.frame(
  Gene = unname(as.character(modules[[gene_col]])),
  Module = unname(as.character(modules[[color_col]])),
  stringsAsFactors = FALSE
)

# 去掉 grey 模块通常更适合下游功能分析；但这里先完整保存。
write.csv(
  gene_module_df,
  file = file.path(FILES_DIR, "Module_Genes_hdWGCNA_signed.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("模块基因列表已保存至：", file.path(FILES_DIR, "Module_Genes_hdWGCNA_signed.csv"), "\n")

# 提取目标模块
target_gene_df <- gene_module_df[gene_module_df$Module %in% TARGET_COLORS, ]

write.csv(
  target_gene_df,
  file = file.path(FILES_DIR, "Target_Selected_Modules_Genes_hdWGCNA_signed.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (color in TARGET_COLORS) {
  single_mod_genes <- gene_module_df$Gene[gene_module_df$Module == color]
  output_path <- file.path(TARGET_MODULES_DIR, paste0("Module_", color, "_hdWGCNA_signed_genes.csv"))

  if (length(single_mod_genes) > 0) {
    write.csv(
      data.frame(Gene = single_mod_genes),
      file = output_path,
      row.names = FALSE,
      quote = FALSE
    )
  } else {
    warning(paste("警告：当前 hdWGCNA 结果中未检测到", color, "模块的基因。"))
  }
}

# hub genes
cat("正在导出 hub genes...\n")
hub_genes <- tryCatch(
  GetHubGenes(malig_seurat, n_hubs = 25, wgcna_name = WGCNA_NAME),
  error = function(e) {
    warning("GetHubGenes 失败：", conditionMessage(e))
    NULL
  }
)

if (!is.null(hub_genes)) {
  write.csv(
    hub_genes,
    file = file.path(FILES_DIR, "Hub_Genes_hdWGCNA_signed_top25.csv"),
    row.names = FALSE,
    quote = FALSE
  )
}

# MEs
MEs <- tryCatch(
  GetMEs(malig_seurat, wgcna_name = WGCNA_NAME),
  error = function(e) {
    warning("GetMEs 失败：", conditionMessage(e))
    NULL
  }
)

if (!is.null(MEs)) {
  write.csv(
    MEs,
    file = file.path(FILES_DIR, "Module_Eigengenes_hdWGCNA_signed.csv"),
    row.names = TRUE,
    quote = FALSE
  )
}

# ==============================================================================
# 12. 保存对象
# ==============================================================================

cat("正在保存 hdWGCNA Seurat 对象与核心结果...\n")

qsave(
  malig_seurat,
  file.path(QS_DIR, "Malignant_seurat_with_hdWGCNA_signed.qs")
)

hdwgcna_results_list <- list(
  wgcna_name = WGCNA_NAME,
  intersect_genes = intersect_genes,
  modules = modules,
  gene_module_df = gene_module_df,
  hub_genes = hub_genes,
  MEs = MEs,
  params = list(
    assay = ASSAY,
    cell_group_col = CELL_GROUP_COL,
    network_group_col = NETWORK_GROUP_COL,
    network_group_name = NETWORK_GROUP_NAME,
    group_by_cols = GROUP_BY_COLS,
    network_type = NETWORK_TYPE,
    tom_type = TOM_TYPE,
    cor_method = COR_METHOD,
    cor_fnc = COR_FNC,
    cor_options = COR_OPTIONS,
    metacell_k = METACELL_K,
    metacell_max_shared = METACELL_MAX_SHARED,
    metacell_min_cells = METACELL_MIN_CELLS,
    metacell_target = METACELL_TARGET
  )
)

qsave(
  hdwgcna_results_list,
  file.path(QS_DIR, "hdWGCNA_Final_Results_signed.qs")
)

cat("全部流程运行完毕！\n")
cat("主要输出：\n")
cat("1. ", file.path(FILES_DIR, "Module_Genes_hdWGCNA_signed.csv"), "\n")
cat("2. ", file.path(FILES_DIR, "Module_Genes_hdWGCNA_signed_full.csv"), "\n")
cat("3. ", file.path(FILES_DIR, "Target_Selected_Modules_Genes_hdWGCNA_signed.csv"), "\n")
cat("4. ", file.path(FILES_DIR, "Hub_Genes_hdWGCNA_signed_top25.csv"), "\n")
cat("5. ", file.path(QS_DIR, "Malignant_seurat_with_hdWGCNA_signed.qs"), "\n")
cat("6. ", file.path(QS_DIR, "hdWGCNA_Final_Results_signed.qs"), "\n")
# ==============================================================================
# 13. 按指定模块颜色提取基因并输出 CSV
# ==============================================================================

cat("正在根据 SELECT_MODULE_COLORS 提取指定模块基因...\n")
# ==============================================================================
# 13 模块颜色提取参数
# ==============================================================================

SELECT_MODULE_COLORS <- c("brown")  # 这里改成你想提取的模块颜色
SELECT_MODULES_DIR <- file.path(FILES_DIR, "Selected_Modules_By_Color")
dir.create(SELECT_MODULES_DIR, showWarnings = FALSE, recursive = TRUE)

selected_module_df <- gene_module_df[
  gene_module_df$Module %in% SELECT_MODULE_COLORS,
]

write.csv(
  selected_module_df,
  file = file.path(FILES_DIR, "Selected_Modules_By_Color_Genes.csv"),
  row.names = FALSE,
  quote = FALSE
)

for (module_color in SELECT_MODULE_COLORS) {
  module_genes <- gene_module_df$Gene[
    gene_module_df$Module == module_color
  ]

  output_file <- file.path(
    SELECT_MODULES_DIR,
    paste0("Module_", module_color, "_genes.csv")
  )

  if (length(module_genes) > 0) {
    write.csv(
      data.frame(Gene = module_genes),
      file = output_file,
      row.names = FALSE,
      quote = FALSE
    )

    cat("已输出模块 ", module_color, " 的基因数：", length(module_genes), "\n")
  } else {
    warning(paste0("未检测到模块颜色：", module_color))
  }
}

cat(
  "指定模块基因已保存至：",
  file.path(FILES_DIR, "Selected_Modules_By_Color_Genes.csv"),
  "\n"
)