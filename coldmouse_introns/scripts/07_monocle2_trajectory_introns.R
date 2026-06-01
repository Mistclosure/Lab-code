# ==============================================================================
# 07_monocle2_trajectory_introns.R
# coldmouse_introns: Monocle 2 拟时序分析 (仅 PBMC Monocyte + Aorta Macrophage)
# ==============================================================================
# 输入：coldmouse_introns_sc_by_tissue_no_harmony_louvain.qs
# 输出：PBMC Monocyte / Aorta Macrophage 的 CDS、伪时间、轨迹图
# 注意：只使用 monocle2，不使用 monocle3
# ==============================================================================

PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))

config <- load_config()
setup_performance(config)

library(Seurat)
library(qs)
library(ggplot2)
library(dplyr)

# ==============================================================================
# monocle2 兼容性补丁 (dplyr / igraph)
# ==============================================================================

# 全局 shim: monocle2 使用了已废弃的 group_by_() 等函数
for (fn_name in c("group_by_", "summarise_", "mutate_", "filter_", "arrange_",
                  "select_", "slice_", "distinct_")) {
  if (!exists(fn_name, envir = globalenv())) {
    base_fn <- gsub("_$", "", fn_name)
    fn_def <- eval(parse(text = paste0("function(.data, ...) { dplyr::", base_fn, "(.data, ...) }")))
    assign(fn_name, fn_def, envir = globalenv())
  }
}

library(monocle)

# 补丁 igraph::dfs 的 neimode 参数问题
if (requireNamespace("igraph", quietly = TRUE)) {
  orig_dfs <- igraph::dfs
  patched_dfs <- function(graph, root, neimode = NULL, mode = NULL, ...) {
    if (!is.null(neimode) && is.null(mode)) mode <- neimode
    orig_dfs(graph, root, mode = mode, ...)
  }
  environment(patched_dfs) <- environment(orig_dfs)
  assignInNamespace("dfs", patched_dfs, ns = "igraph")
}

# 补丁 dplyr 的 deprecated 函数 (monocle2 内部调用)
dplyr_ns <- getNamespace("dplyr")
deprecated_fns <- list(
  select_    = function(.data, ...) select(.data, ...),
  filter_    = function(.data, ...) filter(.data, ...),
  mutate_    = function(.data, ...) mutate(.data, ...),
  arrange_   = function(.data, ...) arrange(.data, ...),
  summarise_ = function(.data, ...) summarise(.data, ...),
  summarize_ = function(.data, ...) summarise(.data, ...),
  group_by_  = function(.data, ...) group_by(.data, ...),
  slice_     = function(.data, ...) slice(.data, ...),
  distinct_  = function(.data, ...) distinct(.data, ...),
  rename_    = function(.data, ...) rename(.data, ...),
  count_     = function(.data, ...) count(.data, ...),
  do_        = function(.data, ...) do(.data, ...),
  pull_      = function(.data, ...) pull(.data, ...),
  transmute_ = function(.data, ...) transmute(.data, ...)
)
for (fn_name in names(deprecated_fns)) {
  tryCatch({
    assignInNamespace(fn_name, deprecated_fns[[fn_name]], ns = "dplyr")
  }, error = function(e) NULL)
}

# ==============================================================================
# 输出目录
# ==============================================================================
output_dir  <- config$paths$output_dir
objects_dir <- file.path(output_dir, "results", "objects")
plots_dir   <- file.path(output_dir, "plots")
files_dir   <- file.path(output_dir, "files")
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(files_dir,   showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 核心函数: run_monocle2_subset_trajectory
# ==============================================================================
run_monocle2_subset_trajectory <- function(obj, subset_name, label_field,
                                           label_value, output_tag) {
  cat("\n===================================================================\n")
  cat("Monocle2 拟时序分析:", output_tag, "\n")
  cat("===================================================================\n")

  # --- 1. 子集检查 ---
  labels <- FetchData(obj, vars = label_field)[[1]]
  if (!label_value %in% unique(labels)) {
    cat("  错误: 在", label_field, "中未找到 '", label_value, "'\n")
    cat("  可用标签:\n")
    print(table(labels))
    stop(label_value, " not found in ", label_field)
  }
  subset_obj <- subset(obj, cells = colnames(obj)[labels == label_value])
  cat("  子集细胞数:", ncol(subset_obj), "\n")
  cat("  Group 分布:\n")
  print(table(subset_obj$Group))
  cat("  seurat_clusters 分布:\n")
  print(table(subset_obj$seurat_clusters))

  # JoinLayers (Seurat v5)
  subset_obj <- JoinLayers(subset_obj)

  # --- 2. 构建 monocle2 CellDataSet ---
  cat("  构建 CellDataSet...\n")
  counts <- GetAssayData(subset_obj, assay = "RNA", layer = "counts")
  pd <- new("AnnotatedDataFrame", data = subset_obj@meta.data)
  fd <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  fd <- new("AnnotatedDataFrame", data = fd)

  cds <- newCellDataSet(
    as(counts, "sparseMatrix"),
    phenoData = pd,
    featureData = fd,
    lowerDetectionLimit = 0.5,
    expressionFamily = negbinomial.size()
  )

  # --- 3. estimateSizeFactors + estimateDispersions ---
  cat("  estimateSizeFactors...\n")
  tryCatch({
    cds <- estimateSizeFactors(cds)
  }, error = function(e) {
    cat("    estimateSizeFactors 失败，手动计算:", e$message, "\n")
    cell_totals <- Matrix::colSums(exprs(cds))
    pData(cds)$Size_Factor <- cell_totals / median(cell_totals)
  })

  cat("  estimateDispersions...\n")
  tryCatch({
    cds <- estimateDispersions(cds)
  }, error = function(e) {
    cat("    estimateDispersions 失败:", e$message, "\n")
    cat("    继续尝试后续步骤...\n")
  })

  # --- 4. Ordering genes: 用 Seurat cluster markers ---
  cat("  计算 Seurat cluster markers 作为 ordering genes...\n")
  Idents(subset_obj) <- "seurat_clusters"
  all_markers <- FindAllMarkers(subset_obj,
                                only.pos = TRUE,
                                min.pct = 0.25,
                                logfc.threshold = 0.25,
                                verbose = FALSE)

  sig_markers <- all_markers %>% filter(p_val_adj < 0.05)
  cat("    显著 marker 数 (p_val_adj < 0.05):", nrow(sig_markers), "\n")

  if (nrow(sig_markers) >= 100) {
    ordering_genes <- unique(sig_markers$gene)
    cat("    使用显著 marker 作为 ordering genes:", length(ordering_genes), "\n")
  } else {
    cat("    显著 marker < 100，放宽为每 cluster top 50...\n")
    top50_per_cluster <- all_markers %>%
      group_by(cluster) %>%
      slice_max(n = 50, order_by = avg_log2FC, with_ties = FALSE)
    ordering_genes <- unique(top50_per_cluster$gene)
    cat("    放宽后 ordering genes:", length(ordering_genes), "\n")

    if (length(ordering_genes) < 100) {
      cat("    仍少于 100 个 ordering genes，停止拟时序分析。\n")
      write.csv(all_markers,
                file.path(files_dir, paste0("introns_monocle2_ordering_genes_", output_tag, ".csv")),
                row.names = FALSE)
      return(invisible(NULL))
    }
  }

  # 过滤掉不在 CDS 中的基因
  ordering_genes <- ordering_genes[ordering_genes %in% rownames(cds)]
  cat("    最终 ordering genes (在 CDS 中):", length(ordering_genes), "\n")

  write.csv(data.frame(gene = ordering_genes, stringsAsFactors = FALSE),
            file.path(files_dir, paste0("introns_monocle2_ordering_genes_", output_tag, ".csv")),
            row.names = FALSE)

  # --- 5. setOrderingFilter + reduceDimension + orderCells ---
  cat("  setOrderingFilter...\n")
  cds <- setOrderingFilter(cds, ordering_genes)

  cat("  reduceDimension (DDRTree)...\n")
  tryCatch({
    cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
  }, error = function(e) {
    cat("    reduceDimension 失败:", e$message, "\n")
    stop("reduceDimension failed for ", output_tag)
  })

  cat("  orderCells...\n")
  tryCatch({
    cds <- orderCells(cds)
  }, error = function(e) {
    cat("    orderCells 失败:", e$message, "\n")
    stop("orderCells failed for ", output_tag)
  })

  # 输出 State 分布统计
  cat("  State x Group:\n")
  print(table(State = pData(cds)$State, Group = pData(cds)$Group))
  cat("  State x seurat_clusters:\n")
  print(table(State = pData(cds)$State, Cluster = pData(cds)$seurat_clusters))

  state_table <- data.frame(
    State = pData(cds)$State,
    Group = pData(cds)$Group,
    seurat_clusters = pData(cds)$seurat_clusters,
    Pseudotime = pData(cds)$Pseudotime
  )
  write.csv(state_table,
            file.path(files_dir, paste0("introns_monocle2_state_table_", output_tag, ".csv")),
            row.names = FALSE)

  # --- 6. 保存 pseudotime metadata ---
  pseudotime_meta <- data.frame(
    Cell = rownames(pData(cds)),
    State = pData(cds)$State,
    Pseudotime = pData(cds)$Pseudotime,
    Group = pData(cds)$Group,
    seurat_clusters = pData(cds)$seurat_clusters,
    stringsAsFactors = FALSE
  )
  write.csv(pseudotime_meta,
            file.path(files_dir, paste0("introns_monocle2_pseudotime_metadata_", output_tag, ".csv")),
            row.names = FALSE)

  # --- 7. 保存 CDS 对象 ---
  cat("  保存 CDS 对象...\n")
  qsave(cds, file.path(objects_dir, paste0("introns_monocle2_cds_", output_tag, ".qs")))
  saveRDS(cds, file.path(objects_dir, paste0("introns_monocle2_cds_", output_tag, ".rds")))

  # --- 8. 可视化 ---
  cat("  生成轨迹图...\n")

  plot_configs <- list(
    list(color_by = "Pseudotime",        w = 8,  h = 6),
    list(color_by = "State",             w = 8,  h = 6),
    list(color_by = "Group",             w = 8,  h = 6),
    list(color_by = "seurat_clusters",   w = 10, h = 6)
  )
  for (pc in plot_configs) {
    tryCatch({
      p <- plot_cell_trajectory(cds, color_by = pc$color_by) +
        ggtitle(paste0("Introns ", output_tag, " Trajectory (", pc$color_by, ")"))
      ggsave(file.path(plots_dir, paste0("introns_monocle2_trajectory_", pc$color_by, "_", output_tag, ".png")),
             plot = p, width = pc$w, height = pc$h, dpi = 300)
      ggsave(file.path(plots_dir, paste0("introns_monocle2_trajectory_", pc$color_by, "_", output_tag, ".pdf")),
             plot = p, width = pc$w, height = pc$h)
    }, error = function(e) {
      cat("    ", pc$color_by, " 图失败:", e$message, "\n")
    })
  }

  cat("  ", output_tag, " 完成\n")
  invisible(cds)
}

# ==============================================================================
# 读取输入对象
# ==============================================================================
sc_by_tissue_file <- file.path(output_dir,
                               paste0("coldmouse_introns_sc_by_tissue_no_harmony_",
                                      config$clustering$method, ".qs"))
cat("读取:", sc_by_tissue_file, "\n")
sc_by_tissue <- qread(sc_by_tissue_file)
cat("可用组织:", paste(names(sc_by_tissue), collapse = ", "), "\n")

# ==============================================================================
# 运行两个子集的拟时序分析
# ==============================================================================

# 1. PBMC Monocyte
run_monocle2_subset_trajectory(
  obj         = sc_by_tissue[["PBMC"]],
  subset_name = "Monocyte",
  label_field = "SingleR.labels",
  label_value = "Monocytes",
  output_tag  = "PBMC_Monocyte"
)

# 2. Aorta Macrophage
run_monocle2_subset_trajectory(
  obj         = sc_by_tissue[["Aorta"]],
  subset_name = "Macrophage",
  label_field = "SingleR.labels",
  label_value = "Macrophages",
  output_tag  = "Aorta_Macrophage"
)

cat("\n=== Monocle2 拟时序分析完成 ===\n")
