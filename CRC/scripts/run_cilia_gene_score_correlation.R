suppressPackageStartupMessages({
  required_pkgs <- c(
    "readxl", "dplyr", "tidyr", "tibble", "ggplot2", "readr",
    "Matrix", "pheatmap", "openxlsx", "Seurat", "qs"
  )
  missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R packages: ", paste(missing_pkgs, collapse = ", "),
      "\nInstall them before running this script."
    )
  }
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(Matrix)
  library(Seurat)
  library(qs)
})

source("/home/zerui/code/pseudo_cilia_score/R/pseudo_cilia_score.R")

# Optional xlsx input. If this file exists, Sheet1/Sheet2 are used.
# Otherwise the project CSV files below are used as the authoritative gene sets.
cilia_xlsx <- NA_character_

gene_files <- c(
  Sheet1 = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table1_primary_cilia_single_gene_confirmed.csv",
  Sheet2 = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table2_primary_cilia_single_gene_confirmed.csv"
)

dataset_dirs <- c(
  GSE132465 = "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE132465",
  GSE178318 = "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE178318",
  GSE225857 = "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE225857",
  GSE231559 = "/mnt/disk1/qiuzerui/downloads/CRC/pseudo_cilia_score_results/GSE231559"
)

outdir <- "/mnt/disk1/qiuzerui/downloads/CRC/cilia_gene_pseudoCiliaScore_correlation_metastasis_direction_4CRCdatasets"
plot_dir <- file.path(outdir, "plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

group_map <- list(
  primary = c("Primary", "primary", "Tumor", "PT", "CRC", "P", "0"),
  metastasis = c("Metastasis", "metastasis", "Met", "LiverMet", "LM", "M", "1")
)

result_dir_by_geneset <- c(
  Sheet1 = "table1_primary_cilia_single_gene_confirmed",
  Sheet2 = "table2_primary_cilia_single_gene_confirmed"
)

gene_col_candidates <- c("gene", "Gene", "symbol", "Symbol", "gene_symbol", "GeneSymbol")

read_gene_column <- function(df) {
  hit <- intersect(gene_col_candidates, colnames(df))
  gene_col <- if (length(hit) > 0) hit[1] else colnames(df)[1]
  genes <- unique(trimws(as.character(df[[gene_col]])))
  genes[!is.na(genes) & genes != ""]
}

read_gene_sets <- function(cilia_xlsx, gene_files) {
  if (!is.na(cilia_xlsx) && file.exists(cilia_xlsx)) {
    sheets <- readxl::excel_sheets(cilia_xlsx)
    use_sheets <- intersect(c("Sheet1", "Sheet2"), sheets)
    if (length(use_sheets) == 0) {
      use_sheets <- sheets[seq_len(min(2, length(sheets)))]
    }
    gene_sets <- lapply(use_sheets, function(sheet) {
      read_gene_column(as.data.frame(readxl::read_excel(cilia_xlsx, sheet = sheet)))
    })
    names(gene_sets) <- paste0("Sheet", seq_along(gene_sets))
    return(gene_sets)
  }

  gene_sets <- lapply(gene_files, function(path) {
    if (!file.exists(path)) {
      stop("Gene set file not found: ", path)
    }
    read_gene_column(read.csv(path, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM"))
  })
  gene_sets
}

standardize_group <- function(x, group_map) {
  x0 <- trimws(as.character(x))
  out <- rep(NA_character_, length(x0))
  out[tolower(x0) %in% tolower(group_map$primary)] <- "Primary"
  out[tolower(x0) %in% tolower(group_map$metastasis)] <- "Metastasis"
  out
}

rename_score_cols <- function(score) {
  pseudo_candidates <- c("pseudo_id", "pseudocell_id", "cell_id", "barcode")
  score_candidates <- c("CiliaScore", "score", "cilia_score")
  pseudo_col <- intersect(pseudo_candidates, colnames(score))[1]
  score_col <- intersect(score_candidates, colnames(score))[1]
  if (is.na(pseudo_col) || is.na(score_col)) {
    stop("score.csv must contain pseudo_id-like and CiliaScore-like columns.")
  }
  if (pseudo_col != "pseudo_id") score$pseudo_id <- score[[pseudo_col]]
  if (score_col != "CiliaScore") score$CiliaScore <- score[[score_col]]
  score
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(list(estimate = NA_real_, p.value = NA_real_, reason = "too_few_or_constant_values"))
  }
  res <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(estimate = NA_real_, p.value = NA_real_, reason = res$message))
  }
  list(estimate = unname(res$estimate), p.value = res$p.value, reason = NA_character_)
}

safe_wilcox <- function(x, group) {
  ok <- is.finite(x) & !is.na(group)
  x <- x[ok]
  group <- group[ok]
  if (length(unique(group)) != 2 || any(table(group) < 1) || length(unique(x)) < 2) {
    return(list(p.value = NA_real_, reason = "insufficient_groups_or_constant_values"))
  }
  res <- tryCatch(
    suppressWarnings(stats::wilcox.test(x ~ group, exact = FALSE)),
    error = function(e) e
  )
  if (inherits(res, "error")) {
    return(list(p.value = NA_real_, reason = res$message))
  }
  list(p.value = res$p.value, reason = NA_character_)
}

load_dataset_object <- function(dataset) {
  if (dataset == "GSE132465") {
    work_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465"
    obj <- qs::qread(file.path(work_dir, "qs", "Seurat", "Malignant_RNA_assay.qs"))
    cli <- read.csv(file.path(work_dir, "metadata", "GSE132465_Cli.csv"), header = TRUE, check.names = FALSE)
    meta <- obj@meta.data
    meta$cell_id_tmp <- rownames(meta)
    meta_merged <- merge(meta, cli, by.x = "orig.ident", by.y = "Tumor", all.x = TRUE, sort = FALSE)
    rownames(meta_merged) <- meta_merged$cell_id_tmp
    meta_merged$cell_id_tmp <- NULL
    obj@meta.data <- meta_merged[colnames(obj), , drop = FALSE]
    tnm_stage <- trimws(obj@meta.data[["TNM stage"]])
    obj$group <- ifelse(grepl("M0$", tnm_stage), "Primary", "Metastasis")
    return(list(obj = obj, sample_col = "orig.ident", group_col = "group"))
  }

  if (dataset == "GSE178318") {
    work_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE178318"
    obj <- qs::qread(file.path(work_dir, "Malignant.qs"))
    obj$type <- sub(".*_", "", colnames(obj))
    obj <- subset(obj, subset = type %in% c("CRC", "LM"))
    obj$group <- ifelse(obj$type == "LM", "Metastasis", "Primary")
    return(list(obj = obj, sample_col = "orig.ident", group_col = "group"))
  }

  if (dataset == "GSE225857") {
    work_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE225857"
    obj <- qs::qread(file.path(work_dir, "Malignant.qs"))
    cli <- read.csv(file.path(work_dir, "GSE225857_Cli.csv"), header = TRUE, check.names = FALSE)
    meta <- obj@meta.data
    meta$cell_id_tmp <- rownames(meta)
    meta_merged <- merge(meta, cli, by.x = "orig.ident", by.y = "Patient ID", all.x = TRUE, sort = FALSE)
    rownames(meta_merged) <- meta_merged$cell_id_tmp
    meta_merged$cell_id_tmp <- NULL
    obj@meta.data <- meta_merged[colnames(obj), , drop = FALSE]
    # Keep the original lowercase labels used when pseudo_id values were created.
    obj$group <- ifelse(obj@meta.data[["Liver metastasis (LM)"]] == "Yes", "metastasis", "primary")
    return(list(obj = obj, sample_col = "orig.ident", group_col = "group"))
  }

  if (dataset == "GSE231559") {
    work_dir <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE231559"
    obj <- qs::qread(file.path(work_dir, "Malignant_RNA_assay.qs"))
    cli <- read.csv(file.path(work_dir, "GSE231559_Cli.csv"), header = TRUE, check.names = FALSE)
    meta <- obj@meta.data
    meta$cell_id_tmp <- rownames(meta)
    meta_merged <- merge(meta, cli, by = "Patient", all.x = TRUE, sort = FALSE)
    rownames(meta_merged) <- meta_merged$cell_id_tmp
    meta_merged$cell_id_tmp <- NULL
    obj@meta.data <- meta_merged[colnames(obj), , drop = FALSE]
    obj$group <- as.character(obj@meta.data[["metastasis"]])
    return(list(obj = obj, sample_col = "orig.ident", group_col = "group"))
  }

  stop("Unsupported dataset: ", dataset)
}

load_score_result <- function(dataset_dir, gene_set_name) {
  result_subdir <- result_dir_by_geneset[[gene_set_name]]
  if (is.null(result_subdir)) {
    stop("No result directory mapping for gene set: ", gene_set_name)
  }
  result_dir <- file.path(dataset_dir, result_subdir)
  score_path <- file.path(result_dir, "score.csv")
  result_path <- file.path(result_dir, "result.rds")
  if (!file.exists(score_path)) stop("Missing score.csv: ", score_path)
  if (!file.exists(result_path)) stop("Missing result.rds: ", result_path)
  list(
    result_dir = result_dir,
    score = rename_score_cols(read.csv(score_path, check.names = FALSE)),
    result = readRDS(result_path)
  )
}

get_pseudo_norm <- function(dataset, gene_set_name, genes, score, result_obj, dataset_cache) {
  if (!is.null(result_obj$pseudo_norm)) {
    pseudo_norm <- result_obj$pseudo_norm
    common <- intersect(score$pseudo_id, colnames(pseudo_norm))
    if (length(common) < max(5, nrow(score) * 0.5)) {
      stop(dataset, " ", gene_set_name, ": too few pseudo_norm columns overlap with score.csv.")
    }
    return(pseudo_norm[, score$pseudo_id[score$pseudo_id %in% common], drop = FALSE])
  }

  ds <- dataset_cache[[dataset]]
  params <- result_obj$parameters
  counts_meta <- extract_counts_and_meta(ds$obj, assay = params$assay %||% "RNA", layer = params$layer %||% "counts")
  pseudo <- build_pseudocells(
    counts = counts_meta$counts,
    meta = counts_meta$meta,
    sample_col = ds$sample_col,
    group_col = ds$group_col,
    celltype_col = NULL,
    subcluster_col = NULL,
    pseudo_mode = params$pseudo_mode %||% "chunk",
    cells_per_pseudo = params$cells_per_pseudo %||% 80,
    min_cells_per_pseudo = params$min_cells_per_pseudo %||% 30,
    seed = params$seed %||% 123,
    verbose = FALSE
  )

  pseudo_counts <- filter_genes_for_scoring(
    pseudo_counts = pseudo$pseudo_counts,
    cilia_genes = genes,
    min_total_count = params$min_total_count %||% 10,
    min_detect_pseudocells = params$min_detect_pseudocells %||% 3,
    min_detect_rate = params$min_detect_rate %||% 0.01,
    remove_mt = params$remove_mt %||% TRUE,
    remove_ribo = params$remove_ribo %||% FALSE
  )
  pseudo_norm <- normalize_pseudocounts(pseudo_counts, scale.factor = params$scale.factor %||% 10000)

  common <- intersect(score$pseudo_id, colnames(pseudo_norm))
  if (length(common) < max(5, nrow(score) * 0.5)) {
    stop(
      dataset, " ", gene_set_name,
      ": rebuilt pseudo_norm has too few overlapping pseudo_id values with score.csv. ",
      "Check grouping columns and seed."
    )
  }
  pseudo_norm[, score$pseudo_id[score$pseudo_id %in% common], drop = FALSE]
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

align_score_norm <- function(score, pseudo_norm, group_map) {
  common <- intersect(score$pseudo_id, colnames(pseudo_norm))
  if (length(common) < max(5, nrow(score) * 0.5)) {
    stop("Too few overlapping pseudo-cells between score and pseudo_norm.")
  }
  score <- score[match(common, score$pseudo_id), , drop = FALSE]
  pseudo_norm <- pseudo_norm[, common, drop = FALSE]
  score$standard_group <- standardize_group(score$group, group_map)
  list(score = score, pseudo_norm = pseudo_norm)
}

calc_gene_score_cor <- function(pseudo_norm, score, genes, dataset, gene_set_name) {
  rows <- lapply(genes, function(gene) {
    if (!gene %in% rownames(pseudo_norm)) {
      return(data.frame(
        dataset = dataset, gene_set_name = gene_set_name, gene = gene,
        status = "missing", mean_expr = NA_real_, detection_rate = NA_real_,
        rho = NA_real_, p_value = NA_real_, pearson_r = NA_real_,
        pearson_p_value = NA_real_, reason = "gene_not_found",
        stringsAsFactors = FALSE
      ))
    }
    expr <- as.numeric(pseudo_norm[gene, ])
    sp <- safe_cor_test(expr, score$CiliaScore, method = "spearman")
    pe <- safe_cor_test(expr, score$CiliaScore, method = "pearson")
    data.frame(
      dataset = dataset, gene_set_name = gene_set_name, gene = gene,
      status = "detected",
      mean_expr = mean(expr, na.rm = TRUE),
      detection_rate = mean(expr > 0, na.rm = TRUE),
      rho = sp$estimate,
      p_value = sp$p.value,
      pearson_r = pe$estimate,
      pearson_p_value = pe$p.value,
      reason = sp$reason,
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out$direction <- ifelse(
    out$status == "detected" & !is.na(out$FDR) & out$FDR < 0.05 & out$rho > 0, "positive",
    ifelse(out$status == "detected" & !is.na(out$FDR) & out$FDR < 0.05 & out$rho < 0, "negative", "not_significant")
  )
  out
}

calc_metastasis_direction <- function(pseudo_norm, score, genes, dataset, gene_set_name) {
  rows <- lapply(genes, function(gene) {
    if (!gene %in% rownames(pseudo_norm)) {
      return(data.frame(
        dataset = dataset, gene_set_name = gene_set_name, gene = gene,
        status = "missing", mean_expr_primary = NA_real_,
        mean_expr_metastasis = NA_real_, median_expr_primary = NA_real_,
        median_expr_metastasis = NA_real_, delta_mean = NA_real_,
        delta_median = NA_real_, p_value = NA_real_, reason = "gene_not_found",
        stringsAsFactors = FALSE
      ))
    }
    expr <- as.numeric(pseudo_norm[gene, ])
    group <- score$standard_group
    primary <- expr[group == "Primary"]
    metastasis <- expr[group == "Metastasis"]
    wt <- safe_wilcox(expr, factor(group, levels = c("Primary", "Metastasis")))
    data.frame(
      dataset = dataset, gene_set_name = gene_set_name, gene = gene,
      status = "detected",
      mean_expr_primary = mean(primary, na.rm = TRUE),
      mean_expr_metastasis = mean(metastasis, na.rm = TRUE),
      median_expr_primary = median(primary, na.rm = TRUE),
      median_expr_metastasis = median(metastasis, na.rm = TRUE),
      delta_mean = mean(metastasis, na.rm = TRUE) - mean(primary, na.rm = TRUE),
      delta_median = median(metastasis, na.rm = TRUE) - median(primary, na.rm = TRUE),
      p_value = wt$p.value,
      reason = wt$reason,
      stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  out$FDR <- p.adjust(out$p_value, method = "BH")
  out$metastasis_direction <- ifelse(
    out$status == "detected" & !is.na(out$FDR) & out$FDR < 0.05 & out$delta_mean > 0, "up_in_metastasis",
    ifelse(out$status == "detected" & !is.na(out$FDR) & out$FDR < 0.05 & out$delta_mean < 0, "down_in_metastasis", "not_significant")
  )
  out
}

calc_sample_level_cor <- function(pseudo_norm, score, genes, dataset, gene_set_name) {
  if (!"sample_id" %in% colnames(score)) {
    stop(dataset, " ", gene_set_name, ": score.csv lacks sample_id.")
  }
  sample_score <- score %>%
    filter(!is.na(sample_id)) %>%
    group_by(sample_id) %>%
    summarise(sample_CiliaScore = mean(CiliaScore, na.rm = TRUE), .groups = "drop")
  if (nrow(sample_score) < 5) {
    warning(dataset, " ", gene_set_name, ": sample-level correlation may be underpowered (n=", nrow(sample_score), ").")
  }

  rows <- lapply(genes, function(gene) {
    if (!gene %in% rownames(pseudo_norm)) {
      return(data.frame(
        dataset = dataset, gene_set_name = gene_set_name, gene = gene,
        status = "missing", n_samples = nrow(sample_score),
        sample_rho = NA_real_, sample_p_value = NA_real_,
        reason = "gene_not_found", stringsAsFactors = FALSE
      ))
    }
    expr_df <- data.frame(
      sample_id = score$sample_id,
      expr = as.numeric(pseudo_norm[gene, ]),
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(sample_id)) %>%
      group_by(sample_id) %>%
      summarise(sample_mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
      inner_join(sample_score, by = "sample_id")
    sp <- safe_cor_test(expr_df$sample_mean_expr, expr_df$sample_CiliaScore, method = "spearman")
    data.frame(
      dataset = dataset, gene_set_name = gene_set_name, gene = gene,
      status = "detected", n_samples = nrow(expr_df),
      sample_rho = sp$estimate, sample_p_value = sp$p.value,
      reason = sp$reason, stringsAsFactors = FALSE
    )
  })
  out <- bind_rows(rows)
  out$sample_FDR <- p.adjust(out$sample_p_value, method = "BH")
  out$sample_direction <- ifelse(
    out$status == "detected" & !is.na(out$sample_FDR) & out$sample_FDR < 0.05 & out$sample_rho > 0, "positive",
    ifelse(out$status == "detected" & !is.na(out$sample_FDR) & out$sample_FDR < 0.05 & out$sample_rho < 0, "negative", "not_significant")
  )
  out
}

summarize_consistency <- function(cor_long, meta_long, gene_set_name) {
  genes <- unique(c(cor_long$gene[cor_long$gene_set_name == gene_set_name], meta_long$gene[meta_long$gene_set_name == gene_set_name]))
  rows <- lapply(genes, function(gene) {
    cdat <- cor_long %>% filter(gene_set_name == !!gene_set_name, gene == !!gene, status == "detected")
    mdat <- meta_long %>% filter(gene_set_name == !!gene_set_name, gene == !!gene, status == "detected")
    n_detected <- length(unique(cdat$dataset))
    n_pos <- sum(cdat$rho > 0, na.rm = TRUE)
    n_neg <- sum(cdat$rho < 0, na.rm = TRUE)
    n_sig_pos <- sum(cdat$direction == "positive", na.rm = TRUE)
    n_sig_neg <- sum(cdat$direction == "negative", na.rm = TRUE)
    n_up <- sum(mdat$delta_mean > 0, na.rm = TRUE)
    n_down <- sum(mdat$delta_mean < 0, na.rm = TRUE)

    corr_label <- if (n_detected < 2) {
      "insufficient"
    } else if (n_pos >= 3 && n_neg < 2) {
      "stable_positive"
    } else if (n_neg >= 3 && n_pos < 2) {
      "stable_negative"
    } else if (n_up >= 3) {
      "metastasis_up_consistent"
    } else if (n_down >= 3) {
      "metastasis_down_consistent"
    } else {
      "mixed"
    }

    has_stable_positive <- n_detected >= 2 && n_pos >= 3 && n_neg < 2
    has_stable_negative <- n_detected >= 2 && n_neg >= 3 && n_pos < 2
    has_meta_up <- n_up >= 3
    has_meta_down <- n_down >= 3
    combined_label <- if (has_stable_positive && has_meta_up) {
      "pro_cilia_score_candidate"
    } else if (has_stable_negative && has_meta_down) {
      "anti_cilia_score_candidate"
    } else if ((has_stable_positive && has_meta_down) || (has_stable_negative && has_meta_up)) {
      "discordant_candidate"
    } else {
      "unstable"
    }

    data.frame(
      gene = gene, gene_set_name = gene_set_name,
      n_datasets_detected = n_detected,
      n_positive_corr = n_pos,
      n_negative_corr = n_neg,
      n_significant_positive_corr = n_sig_pos,
      n_significant_negative_corr = n_sig_neg,
      n_up_in_metastasis = n_up,
      n_down_in_metastasis = n_down,
      mean_rho = mean(cdat$rho, na.rm = TRUE),
      median_rho = median(cdat$rho, na.rm = TRUE),
      mean_delta = mean(mdat$delta_mean, na.rm = TRUE),
      median_delta = median(mdat$delta_mean, na.rm = TRUE),
      consistency_label = corr_label,
      combined_label = combined_label,
      stringsAsFactors = FALSE
    )
  })
  bind_rows(rows)
}

plot_correlation_bar <- function(cor_df, dataset, gene_set_name) {
  x <- cor_df %>% filter(dataset == !!dataset, gene_set_name == !!gene_set_name, status == "detected")
  if (nrow(x) == 0) return(invisible(NULL))
  x <- x %>% arrange(rho) %>% mutate(gene = factor(gene, levels = gene))
  p <- ggplot(x, aes(x = gene, y = rho, fill = direction)) +
    geom_col(width = 0.8) +
    coord_flip() +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "top") +
    labs(title = paste(dataset, gene_set_name, "gene-score Spearman rho"),
         x = NULL, y = "Spearman rho")
  base <- file.path(plot_dir, paste0(dataset, "_", gene_set_name, "_gene_score_correlation_barplot"))
  ggsave(paste0(base, ".png"), p, width = 7, height = max(5, nrow(x) * 0.16), dpi = 300)
  ggsave(paste0(base, ".pdf"), p, width = 7, height = max(5, nrow(x) * 0.16))
  invisible(p)
}

plot_heatmap <- function(long_df, gene_set_name, value_col, file_stub, title) {
  x <- long_df %>% filter(gene_set_name == !!gene_set_name, status == "detected") %>%
    select(gene, dataset, value = all_of(value_col)) %>%
    distinct()
  if (nrow(x) == 0) return(invisible(NULL))
  mat <- x %>%
    pivot_wider(names_from = dataset, values_from = value) %>%
    as.data.frame()
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  mat <- as.matrix(mat)
  mat <- mat[rowSums(!is.na(mat)) > 0, , drop = FALSE]
  if (nrow(mat) == 0) return(invisible(NULL))
  png_path <- file.path(plot_dir, paste0(gene_set_name, "_", file_stub, ".png"))
  pdf_path <- file.path(plot_dir, paste0(gene_set_name, "_", file_stub, ".pdf"))
  grDevices::png(png_path, width = 1800, height = max(1400, nrow(mat) * 35), res = 200)
  pheatmap::pheatmap(
    mat,
    main = title,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  grDevices::dev.off()
  grDevices::pdf(pdf_path, width = 7, height = max(6, nrow(mat) * 0.18))
  pheatmap::pheatmap(
    mat,
    main = title,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
    na_col = "grey90",
    cluster_rows = FALSE,
    cluster_cols = FALSE
  )
  grDevices::dev.off()
  invisible(mat)
}

plot_scatter_examples <- function(example_store, consistency_df) {
  top_genes <- consistency_df %>%
    filter(combined_label == "pro_cilia_score_candidate" | consistency_label == "stable_positive") %>%
    arrange(desc(mean_rho)) %>%
    pull(gene) %>%
    unique() %>%
    head(6)
  if (length(top_genes) == 0) return(invisible(NULL))
  for (gene in top_genes) {
    x <- bind_rows(lapply(example_store, function(df) df[df$gene == gene, , drop = FALSE]))
    if (nrow(x) == 0) next
    p <- ggplot(x, aes(x = expr, y = CiliaScore, color = standard_group)) +
      geom_point(alpha = 0.55, size = 1) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
      facet_wrap(~ dataset, scales = "free") +
      theme_bw() +
      theme(panel.grid = element_blank(), legend.position = "top") +
      labs(title = paste0(gene, " expression vs pseudoCiliaScore"),
           x = "Pseudo-cell log-normalized expression", y = "PseudoCiliaScore")
    base <- file.path(plot_dir, paste0("scatter_", gene, "_expression_vs_CiliaScore"))
    ggsave(paste0(base, ".png"), p, width = 9, height = 6, dpi = 300)
    ggsave(paste0(base, ".pdf"), p, width = 9, height = 6)
  }
}

gene_sets <- read_gene_sets(cilia_xlsx, gene_files)
message("Loaded gene sets: ", paste(names(gene_sets), vapply(gene_sets, length, integer(1)), sep = "=", collapse = ", "))

dataset_cache <- lapply(names(dataset_dirs), load_dataset_object)
names(dataset_cache) <- names(dataset_dirs)

all_cor <- list()
all_meta <- list()
all_sample <- list()
example_store <- list()

for (dataset in names(dataset_dirs)) {
  for (gene_set_name in names(gene_sets)) {
    message("Processing ", dataset, " / ", gene_set_name)
    io <- load_score_result(dataset_dirs[[dataset]], gene_set_name)
    genes <- gene_sets[[gene_set_name]]
    pseudo_norm <- get_pseudo_norm(dataset, gene_set_name, genes, io$score, io$result, dataset_cache)
    aligned <- align_score_norm(io$score, pseudo_norm, group_map)
    score <- aligned$score
    pseudo_norm <- aligned$pseudo_norm

    cor_df <- calc_gene_score_cor(pseudo_norm, score, genes, dataset, gene_set_name)
    meta_df <- calc_metastasis_direction(pseudo_norm, score, genes, dataset, gene_set_name)
    sample_df <- calc_sample_level_cor(pseudo_norm, score, genes, dataset, gene_set_name)

    readr::write_csv(cor_df, file.path(outdir, paste0(dataset, "_", gene_set_name, "_gene_score_correlation.csv")))
    readr::write_csv(meta_df, file.path(outdir, paste0(dataset, "_", gene_set_name, "_gene_metastasis_direction.csv")))
    readr::write_csv(sample_df, file.path(outdir, paste0(dataset, "_", gene_set_name, "_sample_level_correlation.csv")))

    all_cor[[paste(dataset, gene_set_name, sep = "__")]] <- cor_df
    all_meta[[paste(dataset, gene_set_name, sep = "__")]] <- meta_df
    all_sample[[paste(dataset, gene_set_name, sep = "__")]] <- sample_df

    detected_genes <- genes[genes %in% rownames(pseudo_norm)]
    keep_examples <- head(detected_genes, 10)
    if (length(keep_examples) > 0) {
      example_store[[paste(dataset, gene_set_name, sep = "__")]] <- bind_rows(lapply(keep_examples, function(gene) {
        data.frame(
          dataset = dataset, gene_set_name = gene_set_name, gene = gene,
          expr = as.numeric(pseudo_norm[gene, ]),
          CiliaScore = score$CiliaScore,
          standard_group = score$standard_group,
          stringsAsFactors = FALSE
        )
      }))
    }
  }
}

cor_long <- bind_rows(all_cor)
meta_long <- bind_rows(all_meta)
sample_long <- bind_rows(all_sample)

readr::write_csv(cor_long, file.path(outdir, "all_datasets_all_genes_correlation_long.csv"))
readr::write_csv(meta_long, file.path(outdir, "all_datasets_all_genes_metastasis_direction_long.csv"))
readr::write_csv(sample_long, file.path(outdir, "all_datasets_all_genes_sample_level_correlation_long.csv"))

consistency_list <- lapply(names(gene_sets), function(gs) summarize_consistency(cor_long, meta_long, gs))
names(consistency_list) <- names(gene_sets)
consistency_all <- bind_rows(consistency_list)

for (gs in names(consistency_list)) {
  readr::write_csv(consistency_list[[gs]], file.path(outdir, paste0(gs, "_cross_dataset_gene_consistency.csv")))
}

stable_positive_list <- lapply(names(consistency_list), function(gs) {
  x <- consistency_list[[gs]] %>%
    filter(consistency_label == "stable_positive") %>%
    arrange(desc(mean_rho), gene)
  readr::write_csv(x, file.path(outdir, paste0(gs, "_stable_positive_genes.csv")))
  writeLines(x$gene, file.path(outdir, paste0(gs, "_stable_positive_gene_symbols.txt")))
  x
})
names(stable_positive_list) <- names(consistency_list)
stable_positive_all <- bind_rows(stable_positive_list)
readr::write_csv(stable_positive_all, file.path(outdir, "all_stable_positive_genes.csv"))
writeLines(unique(stable_positive_all$gene), file.path(outdir, "all_stable_positive_gene_symbols.txt"))

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "all_consistency")
openxlsx::writeData(wb, "all_consistency", consistency_all)
for (gs in names(consistency_list)) {
  openxlsx::addWorksheet(wb, gs)
  openxlsx::writeData(wb, gs, consistency_list[[gs]])
}
openxlsx::addWorksheet(wb, "correlation_long")
openxlsx::writeData(wb, "correlation_long", cor_long)
openxlsx::addWorksheet(wb, "metastasis_direction")
openxlsx::writeData(wb, "metastasis_direction", meta_long)
openxlsx::addWorksheet(wb, "stable_positive")
openxlsx::writeData(wb, "stable_positive", stable_positive_all)
for (gs in names(stable_positive_list)) {
  sheet_name <- paste0(gs, "_stable_positive")
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, stable_positive_list[[gs]])
}
openxlsx::saveWorkbook(wb, file.path(outdir, "all_gene_consistency_summary.xlsx"), overwrite = TRUE)

for (dataset in names(dataset_dirs)) {
  for (gene_set_name in names(gene_sets)) {
    plot_correlation_bar(cor_long, dataset, gene_set_name)
  }
}

for (gene_set_name in names(gene_sets)) {
  plot_heatmap(cor_long, gene_set_name, "rho", "correlation_rho_heatmap", paste0(gene_set_name, " Spearman rho"))
  plot_heatmap(meta_long, gene_set_name, "delta_mean", "metastasis_delta_mean_heatmap", paste0(gene_set_name, " Metastasis - Primary delta"))
}

summary_plot_df <- consistency_all %>%
  count(gene_set_name, consistency_label, name = "n_genes")
p_summary <- ggplot(summary_plot_df, aes(x = consistency_label, y = n_genes, fill = consistency_label)) +
  geom_col(width = 0.75) +
  facet_wrap(~ gene_set_name) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Cilia gene consistency summary", x = NULL, y = "Gene count")
ggsave(file.path(plot_dir, "consistency_summary_plot.png"), p_summary, width = 8, height = 5, dpi = 300)
ggsave(file.path(plot_dir, "consistency_summary_plot.pdf"), p_summary, width = 8, height = 5)

plot_scatter_examples(example_store, consistency_all)

message("Finished cilia gene score correlation analysis.")
message("Outputs written to: ", outdir)
