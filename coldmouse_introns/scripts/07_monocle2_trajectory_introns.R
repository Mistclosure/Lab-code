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

# 补丁 monocle2 内部对 igraph::nei() 的旧调用；igraph >= 2.1 已改为 .nei()。
patch_monocle_nei <- function() {
  ns <- getNamespace("monocle")
  for (fn_name in ls(ns)) {
    fn <- get(fn_name, envir = ns)
    if (!is.function(fn)) next
    body_txt <- paste(deparse(body(fn)), collapse = "\n")
    if (!grepl("nei(", body_txt, fixed = TRUE)) next

    patched_txt <- gsub("(?<![.])nei\\(", ".nei(", body_txt, perl = TRUE)
    patched_body <- tryCatch(parse(text = patched_txt)[[1]], error = function(e) NULL)
    if (is.null(patched_body)) next

    body(fn) <- patched_body
    tryCatch({
      assignInNamespace(fn_name, fn, ns = "monocle")
    }, error = function(e) NULL)
  }
}
patch_monocle_nei()

# 补丁 dplyr 的 deprecated 函数 (monocle2 内部调用)
# monocle2 调用了 dplyr 0.7 之前的 *_() SE 接口；新版 dplyr 中这些函数已 defunct。
make_compat_select_ <- function() {
  function(.data, ..., .dots = NULL) {
    lazies <- lazyeval::lazy_dots(...)
    vals <- lapply(lazies, function(z) {
      tryCatch(lazyeval::lazy_eval(z), error = function(e) as.character(z$expr))
    })
    if (!is.null(.dots) && length(.dots) > 0) vals <- c(vals, .dots)
    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0) return(.data)

    cols <- unlist(lapply(vals, function(x) {
      if (length(x) == 0) return(character(0))
      x1 <- x[[1]]
      if (is.numeric(x1)) {
        return(names(.data)[as.integer(x1)])
      }
      as.character(x1)
    }), use.names = FALSE)
    cols <- cols[nzchar(cols)]
    out <- dplyr::select(.data, dplyr::all_of(cols))

    new_names <- names(vals)
    if (!is.null(new_names)) {
      named <- !is.na(new_names) & nzchar(new_names)
      if (any(named)) {
        idx <- match(cols[named], names(out))
        ok <- !is.na(idx)
        names(out)[idx[ok]] <- new_names[named][ok]
      }
    }
    out
  }
}
make_compat_group_by_ <- function() {
  function(.data, ..., .dots = NULL) {
    lazies <- lazyeval::lazy_dots(...)
    vals <- lapply(lazies, function(z) {
      tryCatch(lazyeval::lazy_eval(z), error = function(e) as.character(z$expr))
    })
    if (!is.null(.dots) && length(.dots) > 0) vals <- c(vals, .dots)
    vals <- vals[!vapply(vals, is.null, logical(1))]
    if (length(vals) == 0) return(dplyr::group_by(.data))
    cols <- unlist(lapply(vals, function(x) as.character(x[[1]])), use.names = FALSE)
    cols <- cols[nzchar(cols)]
    dplyr::group_by(.data, dplyr::across(dplyr::all_of(cols)))
  }
}
make_compat_filter_ <- function() {
  function(.data, ..., .dots = NULL) dplyr::filter(.data, ...)
}
make_compat_arrange_ <- function() {
  function(.data, ..., .dots = NULL) dplyr::arrange(.data, ...)
}
make_compat_mutate_ <- function() {
  function(.data, ..., .dots = NULL) dplyr::mutate(.data, ...)
}
make_compat_summarise_ <- function() {
  function(.data, ..., .dots = NULL) dplyr::summarise(.data, ...)
}

patch_binding <- function(env, name, value) {
  if (!exists(name, envir = env, inherits = FALSE)) return(invisible(FALSE))
  was_locked <- bindingIsLocked(name, env)
  if (was_locked) unlockBinding(name, env)
  assign(name, value, envir = env)
  if (was_locked) lockBinding(name, env)
  invisible(TRUE)
}

dplyr_compat_fns <- list(
  select_    = make_compat_select_(),
  group_by_  = make_compat_group_by_(),
  filter_    = make_compat_filter_(),
  arrange_   = make_compat_arrange_(),
  mutate_    = make_compat_mutate_(),
  summarise_ = make_compat_summarise_(),
  summarize_ = make_compat_summarise_()
)

dplyr_ns <- getNamespace("dplyr")
for (fn_name in names(dplyr_compat_fns)) {
  patch_binding(dplyr_ns, fn_name, dplyr_compat_fns[[fn_name]])
  if ("package:dplyr" %in% search()) {
    patch_binding(as.environment("package:dplyr"), fn_name, dplyr_compat_fns[[fn_name]])
  }
  monocle_imports <- parent.env(getNamespace("monocle"))
  patch_binding(monocle_imports, fn_name, dplyr_compat_fns[[fn_name]])
  assign(fn_name, dplyr_compat_fns[[fn_name]], envir = globalenv())
}


# ==============================================================================
# 输出目录
# ==============================================================================
output_dir  <- config$paths$output_dir
objects_dir <- get_results_dir(config, "objects")
rds_dir     <- get_results_dir(config, "rds")
plots_dir   <- get_results_dir(config, "plots")
files_dir   <- get_results_dir(config, "files")
dir.create(objects_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rds_dir,     showWarnings = FALSE, recursive = TRUE)
dir.create(plots_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(files_dir,   showWarnings = FALSE, recursive = TRUE)

# monocle2 DDRTree 对全量大子集很慢；这里限制进入轨迹建模的规模，
# 保留 Group x seurat_clusters 的组成，避免单一组或 cluster 主导。
MAX_MONOCLE2_CELLS <- 3000
MAX_ORDERING_GENES <- 800

get_monocle2_output_dirs <- function(output_tag) {
  rel_dir <- switch(output_tag,
    PBMC_Monocyte = file.path("PBMC", "Monocytes"),
    Aorta_Macrophage = file.path("Aorta", "Macrophage"),
    output_tag
  )
  tag_files_dir <- file.path(files_dir, rel_dir)
  tag_plots_dir <- file.path(plots_dir, rel_dir)
  dir.create(tag_files_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(tag_plots_dir, showWarnings = FALSE, recursive = TRUE)
  list(files = tag_files_dir, plots = tag_plots_dir)
}

set.seed(20260602)

stratified_downsample_cells <- function(obj, max_cells = MAX_MONOCLE2_CELLS) {
  if (ncol(obj) <= max_cells) {
    return(colnames(obj))
  }

  meta <- obj@meta.data
  meta$.cell <- rownames(meta)
  meta$.stratum <- paste(meta$Group, meta$seurat_clusters, sep = "__")

  strata <- split(meta$.cell, meta$.stratum)
  strata <- strata[lengths(strata) > 0]

  target <- pmax(1, floor(lengths(strata) / ncol(obj) * max_cells))
  overflow <- sum(target) - max_cells
  if (overflow > 0) {
    reducible <- which(target > 1)
    while (overflow > 0 && length(reducible) > 0) {
      i <- reducible[which.max(target[reducible])]
      target[i] <- target[i] - 1
      overflow <- overflow - 1
      reducible <- which(target > 1)
    }
  }

  sampled <- unlist(Map(function(cells, n) {
    sample(cells, min(length(cells), n))
  }, strata, target), use.names = FALSE)

  if (length(sampled) > max_cells) {
    sampled <- sample(sampled, max_cells)
  }
  sampled
}

# ==============================================================================
# 核心函数: run_monocle2_subset_trajectory
# ==============================================================================
run_monocle2_subset_trajectory <- function(obj, subset_name, label_field,
                                           label_value, output_tag) {
  out_dirs <- get_monocle2_output_dirs(output_tag)
  tag_files_dir <- out_dirs$files
  tag_plots_dir <- out_dirs$plots

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
  cat("  原始子集细胞数:", ncol(subset_obj), "\n")
  cat("  Group 分布:\n")
  print(table(subset_obj$Group))
  cat("  seurat_clusters 分布:\n")
  print(table(subset_obj$seurat_clusters))

  sampled_cells <- stratified_downsample_cells(subset_obj, MAX_MONOCLE2_CELLS)
  if (length(sampled_cells) < ncol(subset_obj)) {
    cat("  monocle2 DDRTree 输入过大，按 Group x seurat_clusters 分层抽样:",
        ncol(subset_obj), "->", length(sampled_cells), "cells\n")
    subset_obj <- subset(subset_obj, cells = sampled_cells)
    cat("  抽样后 Group 分布:\n")
    print(table(subset_obj$Group))
    cat("  抽样后 seurat_clusters 分布:\n")
    print(table(subset_obj$seurat_clusters))
  }

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
    ordering_genes <- sig_markers %>%
      arrange(p_val_adj, desc(avg_log2FC)) %>%
      pull(gene) %>%
      unique()
    if (length(ordering_genes) > MAX_ORDERING_GENES) {
      ordering_genes <- ordering_genes[seq_len(MAX_ORDERING_GENES)]
    }
    cat("    使用显著 marker 作为 ordering genes:", length(ordering_genes), "\n")
  } else {
    cat("    显著 marker < 100，放宽为每 cluster top 50...\n")
    top50_per_cluster <- all_markers %>%
      group_by(cluster) %>%
      slice_max(n = 50, order_by = avg_log2FC, with_ties = FALSE)
    ordering_genes <- top50_per_cluster %>%
      arrange(desc(avg_log2FC)) %>%
      pull(gene) %>%
      unique()
    if (length(ordering_genes) > MAX_ORDERING_GENES) {
      ordering_genes <- ordering_genes[seq_len(MAX_ORDERING_GENES)]
    }
    cat("    放宽后 ordering genes:", length(ordering_genes), "\n")

    if (length(ordering_genes) < 100) {
      cat("    仍少于 100 个 ordering genes，停止拟时序分析。\n")
      write.csv(all_markers,
                file.path(tag_files_dir, paste0("introns_monocle2_ordering_genes_", output_tag, ".csv")),
                row.names = FALSE)
      return(invisible(NULL))
    }
  }

  # 过滤掉不在 CDS 中的基因
  ordering_genes <- ordering_genes[ordering_genes %in% rownames(cds)]
  cat("    最终 ordering genes (在 CDS 中):", length(ordering_genes), "\n")

  write.csv(data.frame(gene = ordering_genes, stringsAsFactors = FALSE),
            file.path(tag_files_dir, paste0("introns_monocle2_ordering_genes_", output_tag, ".csv")),
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
            file.path(tag_files_dir, paste0("introns_monocle2_state_table_", output_tag, ".csv")),
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
            file.path(tag_files_dir, paste0("introns_monocle2_pseudotime_metadata_", output_tag, ".csv")),
            row.names = FALSE)

  # --- 7. 保存 CDS 对象 ---
  cat("  保存 CDS 对象...\n")
  qsave(cds, file.path(objects_dir, paste0("introns_monocle2_cds_", output_tag, ".qs")))
  saveRDS(cds, file.path(rds_dir, paste0("introns_monocle2_cds_", output_tag, ".rds")))

  # --- 8. 可视化：使用 Seurat UMAP 坐标展示 monocle2 结果 ---
  cat("  生成 UMAP 拟时序图...\n")

  if (!"umap" %in% names(subset_obj@reductions)) {
    stop("Seurat subset has no UMAP reduction for ", output_tag)
  }
  umap_coords <- Embeddings(subset_obj, "umap")
  cds_cells <- rownames(pData(cds))
  missing_umap_cells <- setdiff(cds_cells, rownames(umap_coords))
  if (length(missing_umap_cells) > 0) {
    stop("UMAP coordinates missing for ", length(missing_umap_cells), " cells in ", output_tag)
  }

  umap_df <- data.frame(
    Cell = cds_cells,
    UMAP_1 = umap_coords[cds_cells, 1],
    UMAP_2 = umap_coords[cds_cells, 2],
    Pseudotime = pData(cds)$Pseudotime,
    State = pData(cds)$State,
    Group = pData(cds)$Group,
    seurat_clusters = pData(cds)$seurat_clusters,
    stringsAsFactors = FALSE
  )
  write.csv(umap_df,
            file.path(tag_files_dir, paste0("introns_monocle2_umap_metadata_", output_tag, ".csv")),
            row.names = FALSE)

  plot_monocle2_umap <- function(df, color_by, output_tag) {
    plot_df <- df
    base_plot <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]])) +
      geom_point(size = 0.35, alpha = 0.85) +
      labs(
        title = paste0("Introns ", output_tag, " monocle2 on UMAP (", color_by, ")"),
        x = "UMAP_1",
        y = "UMAP_2",
        color = color_by
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8)
      )

    if (is.numeric(plot_df[[color_by]])) {
      base_plot + scale_color_gradientn(
        colors = c("#2c7bb6", "#ffffbf", "#d7191c"),
        na.value = "grey80"
      )
    } else {
      plot_df[[color_by]] <- factor(plot_df[[color_by]])
      ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]])) +
        geom_point(size = 0.35, alpha = 0.85) +
        labs(
          title = paste0("Introns ", output_tag, " monocle2 on UMAP (", color_by, ")"),
          x = "UMAP_1",
          y = "UMAP_2",
          color = color_by
        ) +
        theme_classic(base_size = 12) +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 8)
        )
    }
  }

  plot_configs <- list(
    list(color_by = "Pseudotime",        w = 8,  h = 6),
    list(color_by = "State",             w = 8,  h = 6),
    list(color_by = "Group",             w = 8,  h = 6),
    list(color_by = "seurat_clusters",   w = 10, h = 6)
  )
  for (pc in plot_configs) {
    tryCatch({
      p <- plot_monocle2_umap(umap_df, pc$color_by, output_tag)
      ggsave(file.path(tag_plots_dir, paste0("introns_monocle2_umap_", pc$color_by, "_", output_tag, ".png")),
             plot = p, width = pc$w, height = pc$h, dpi = 300)
      ggsave(file.path(tag_plots_dir, paste0("introns_monocle2_umap_", pc$color_by, "_", output_tag, ".pdf")),
             plot = p, width = pc$w, height = pc$h)
    }, error = function(e) {
      cat("    ", pc$color_by, " UMAP 图失败:", e$message, "\n")
    })
  }

  cat("  ", output_tag, " 完成\n")
  invisible(cds)
}

# ==============================================================================
# 读取输入对象
# ==============================================================================
sc_by_tissue_file <- file.path(objects_dir,
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
