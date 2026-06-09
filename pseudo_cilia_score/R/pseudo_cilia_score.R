# Pseudo-cell-level cilia activity scoring
# Adapted from background-corrected module scoring strategies.

extract_counts_and_meta <- function(object,
                                    meta = NULL,
                                    assay = "RNA",
                                    layer = "counts") {
  is_seurat <- inherits(object, "Seurat")

  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("SeuratObject is required to extract counts from a Seurat object.")
    }
    counts <- tryCatch(
      SeuratObject::LayerData(object = object, assay = assay, layer = layer),
      error = function(e) {
        SeuratObject::GetAssayData(object = object, assay = assay, slot = layer)
      }
    )
    meta <- object@meta.data
  } else {
    if (!(is.matrix(object) || inherits(object, "Matrix"))) {
      stop("object must be a Seurat object, matrix, or Matrix sparse matrix.")
    }
    counts <- object
    if (is.null(meta)) {
      stop("meta must be provided when object is a matrix or sparse matrix.")
    }
  }

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have gene rownames and cell colnames.")
  }
  if (is.null(rownames(meta))) {
    stop("metadata must have rownames matching cell names.")
  }

  common_cells <- intersect(colnames(counts), rownames(meta))
  if (length(common_cells) == 0) {
    stop("No overlapping cells between counts colnames and metadata rownames.")
  }
  if (length(common_cells) < ncol(counts) || length(common_cells) < nrow(meta)) {
    warning("Counts and metadata cell names are not identical; using overlapping cells only.")
  }

  counts <- counts[, common_cells, drop = FALSE]
  meta <- meta[common_cells, , drop = FALSE]

  list(counts = counts, meta = meta)
}

subset_target_cells <- function(counts,
                                meta,
                                celltype_col = NULL,
                                target_celltype = NULL) {
  if (is.null(target_celltype)) {
    return(list(counts = counts, meta = meta))
  }
  if (is.null(celltype_col) || !celltype_col %in% colnames(meta)) {
    stop("celltype_col must be provided and present in metadata when target_celltype is set.")
  }

  keep <- !is.na(meta[[celltype_col]]) & meta[[celltype_col]] %in% target_celltype
  if (!any(keep)) {
    stop("No cells remain after target_celltype filtering.")
  }

  list(counts = counts[, keep, drop = FALSE], meta = meta[keep, , drop = FALSE])
}

.make_group_key <- function(meta, grouping_cols) {
  parts <- lapply(grouping_cols, function(x) {
    value <- as.character(meta[[x]])
    value[is.na(value) | value == ""] <- "NA"
    value
  })
  do.call(paste, c(parts, sep = "||"))
}

.make_safe_id <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

build_pseudocells <- function(counts,
                              meta,
                              sample_col = "sample_id",
                              group_col = "group",
                              celltype_col = NULL,
                              subcluster_col = NULL,
                              pseudo_mode = c("chunk", "sample"),
                              cells_per_pseudo = 80,
                              min_cells_per_pseudo = 30,
                              seed = 123,
                              verbose = TRUE) {
  pseudo_mode <- match.arg(pseudo_mode)
  required_cols <- c(sample_col, group_col, celltype_col, subcluster_col)
  required_cols <- required_cols[!is.null(required_cols)]
  missing_cols <- setdiff(required_cols, colnames(meta))
  if (length(missing_cols) > 0) {
    stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))
  }

  grouping_cols <- required_cols
  group_key <- .make_group_key(meta, grouping_cols)
  strata <- split(rownames(meta), group_key)
  set.seed(seed)

  mapping <- list()
  pseudo_meta_list <- list()
  pseudo_index <- 0L

  for (key in names(strata)) {
    cells <- strata[[key]]
    n_cells <- length(cells)
    if (n_cells < min_cells_per_pseudo) {
      warning("Dropping stratum with too few cells: ", key, " (n=", n_cells, ").")
      next
    }

    if (pseudo_mode == "sample") {
      chunks <- list(cells)
    } else {
      shuffled <- sample(cells, length(cells))
      n_chunks <- floor(length(shuffled) / cells_per_pseudo)
      remainder <- length(shuffled) %% cells_per_pseudo

      chunks <- list()
      if (n_chunks > 0) {
        chunks <- split(shuffled[seq_len(n_chunks * cells_per_pseudo)],
                        rep(seq_len(n_chunks), each = cells_per_pseudo))
      }
      if (remainder >= min_cells_per_pseudo) {
        chunks[[length(chunks) + 1L]] <- tail(shuffled, remainder)
      }
      if (length(chunks) == 0L && n_cells >= min_cells_per_pseudo) {
        chunks <- list(shuffled)
      }
    }

    stratum_meta <- meta[cells[1], , drop = FALSE]
    for (i in seq_along(chunks)) {
      pseudo_index <- pseudo_index + 1L
      pseudo_id <- paste0("pseudo_", pseudo_index, "_", .make_safe_id(key), "_", i)
      chunk_cells <- chunks[[i]]
      mapping[[pseudo_id]] <- data.frame(
        cell = chunk_cells,
        pseudo_id = pseudo_id,
        stringsAsFactors = FALSE
      )
      pseudo_meta_list[[pseudo_id]] <- data.frame(
        pseudo_id = pseudo_id,
        sample_id = as.character(stratum_meta[[sample_col]]),
        group = as.character(stratum_meta[[group_col]]),
        cell_type = if (!is.null(celltype_col)) as.character(stratum_meta[[celltype_col]]) else NA_character_,
        subcluster = if (!is.null(subcluster_col)) as.character(stratum_meta[[subcluster_col]]) else NA_character_,
        n_cells = length(chunk_cells),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(mapping) == 0L) {
    stop("No pseudo-cells were constructed. Check min_cells_per_pseudo and grouping columns.")
  }

  cell_map <- do.call(rbind, mapping)
  pseudo_meta <- do.call(rbind, pseudo_meta_list)
  rownames(pseudo_meta) <- pseudo_meta$pseudo_id

  pseudo_levels <- pseudo_meta$pseudo_id
  cell_index <- match(cell_map$cell, colnames(counts))
  pseudo_index <- match(cell_map$pseudo_id, pseudo_levels)

  design <- Matrix::sparseMatrix(
    i = cell_index,
    j = pseudo_index,
    x = 1,
    dims = c(ncol(counts), length(pseudo_levels)),
    dimnames = list(colnames(counts), pseudo_levels)
  )
  pseudo_counts <- counts %*% design
  pseudo_counts <- as(pseudo_counts, "dgCMatrix")

  if (verbose) {
    message("Constructed ", ncol(pseudo_counts), " pseudo-cells from ",
            length(unique(cell_map$cell)), " cells.")
  }

  list(pseudo_counts = pseudo_counts, pseudo_meta = pseudo_meta, cell_map = cell_map)
}

filter_genes_for_scoring <- function(pseudo_counts,
                                     cilia_genes,
                                     min_total_count = 10,
                                     min_detect_pseudocells = 3,
                                     min_detect_rate = 0.01,
                                     remove_mt = TRUE,
                                     remove_ribo = FALSE) {
  total_count <- Matrix::rowSums(pseudo_counts)
  detected_n <- Matrix::rowSums(pseudo_counts > 0)
  detection_rate <- detected_n / ncol(pseudo_counts)

  basic_keep <- total_count >= min_total_count &
    detected_n >= min_detect_pseudocells &
    detection_rate >= min_detect_rate

  unwanted <- rep(FALSE, nrow(pseudo_counts))
  genes <- rownames(pseudo_counts)
  if (remove_mt) {
    unwanted <- unwanted | grepl("^MT-|^mt-", genes)
  }
  if (remove_ribo) {
    unwanted <- unwanted | grepl("^RPL|^RPS|^Rpl|^Rps", genes)
  }

  is_cilia <- genes %in% cilia_genes
  keep <- basic_keep & (!unwanted | is_cilia)
  pseudo_counts[keep, , drop = FALSE]
}

normalize_pseudocounts <- function(pseudo_counts, scale.factor = 10000) {
  lib_size <- Matrix::colSums(pseudo_counts)
  if (any(lib_size <= 0)) {
    stop("At least one pseudo-cell has zero library size after filtering.")
  }
  norm <- t(t(pseudo_counts) / lib_size * scale.factor)
  log1p(norm)
}

.make_bins <- function(values, n_bins) {
  n <- length(values)
  n_bins <- max(1L, min(as.integer(n_bins), n))
  ord <- order(values, seq_along(values), na.last = TRUE)
  bins <- integer(n)
  bins[ord] <- ceiling(seq_along(ord) / n * n_bins)
  bins[bins < 1L] <- 1L
  bins[bins > n_bins] <- n_bins
  bins
}

compute_gene_background_features <- function(pseudo_counts,
                                             pseudo_norm,
                                             cilia_genes_used,
                                             n.mean.bins = 20,
                                             n.detect.bins = 10) {
  n_genes <- nrow(pseudo_norm)
  n.mean.bins <- max(1L, min(as.integer(n.mean.bins), n_genes))
  n.detect.bins <- max(1L, min(as.integer(n.detect.bins), n_genes))

  mean_expr <- Matrix::rowMeans(pseudo_norm)
  detection_rate <- Matrix::rowSums(pseudo_counts > 0) / ncol(pseudo_counts)

  data.frame(
    gene = rownames(pseudo_norm),
    mean_expr = as.numeric(mean_expr),
    detection_rate = as.numeric(detection_rate),
    mean_bin = .make_bins(mean_expr, n.mean.bins),
    detect_bin = .make_bins(detection_rate, n.detect.bins),
    is_cilia_gene = rownames(pseudo_norm) %in% cilia_genes_used,
    used_as_cilia_gene = rownames(pseudo_norm) %in% cilia_genes_used,
    stringsAsFactors = FALSE
  )
}

match_background_genes <- function(gene_info,
                                   cilia_genes_used,
                                   n.control = 50,
                                   verbose = TRUE) {
  candidate_info <- gene_info[!gene_info$gene %in% cilia_genes_used, , drop = FALSE]
  max_mean_radius <- max(gene_info$mean_bin, na.rm = TRUE)
  max_detect_radius <- max(gene_info$detect_bin, na.rm = TRUE)
  max_radius <- max(max_mean_radius, max_detect_radius)

  mapping <- vector("list", length(cilia_genes_used))
  names(mapping) <- cilia_genes_used

  for (gene in cilia_genes_used) {
    row <- gene_info[gene_info$gene == gene, , drop = FALSE]
    if (nrow(row) != 1L) {
      warning("Skipping cilia gene without gene_info row: ", gene)
      next
    }

    selected <- character(0)
    radius_used <- NA_integer_
    for (radius in 0:max_radius) {
      in_bin <- abs(candidate_info$mean_bin - row$mean_bin) <= radius &
        abs(candidate_info$detect_bin - row$detect_bin) <= radius
      selected <- candidate_info$gene[in_bin]
      radius_used <- radius
      if (length(selected) >= n.control) {
        break
      }
    }

    if (length(selected) == 0L) {
      warning("No matched background genes found for ", gene, ".")
    } else if (length(selected) < n.control && verbose) {
      warning("Only ", length(selected), " matched background genes found for ",
              gene, "; using available genes.")
    }

    mapping[[gene]] <- data.frame(
      cilia_gene = gene,
      background_gene = selected,
      mean_bin = row$mean_bin,
      detect_bin = row$detect_bin,
      expansion_radius = radius_used,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, mapping)
}

calculate_score_with_random_controls <- function(pseudo_norm,
                                                 cilia_genes_used,
                                                 background_mapping,
                                                 n.control = 50,
                                                 num.rounds = 1000,
                                                 seed = 123,
                                                 verbose = TRUE) {
  if (length(cilia_genes_used) < 1L) {
    stop("No cilia genes available for scoring.")
  }
  target_score <- Matrix::colMeans(pseudo_norm[cilia_genes_used, , drop = FALSE])

  split_bg <- split(background_mapping$background_gene, background_mapping$cilia_gene)
  split_bg <- split_bg[cilia_genes_used]

  set.seed(seed)
  background_scores <- matrix(NA_real_, nrow = ncol(pseudo_norm), ncol = num.rounds)
  rownames(background_scores) <- colnames(pseudo_norm)

  for (round_i in seq_len(num.rounds)) {
    sampled_genes <- unlist(lapply(split_bg, function(pool) {
      pool <- unique(pool)
      pool <- pool[!is.na(pool) & pool %in% rownames(pseudo_norm)]
      if (length(pool) == 0L) {
        return(character(0))
      }
      sample(pool, size = min(n.control, length(pool)), replace = FALSE)
    }), use.names = FALSE)

    sampled_genes <- unique(sampled_genes)
    if (length(sampled_genes) == 0L) {
      stop("No background genes sampled; check background matching.")
    }
    background_scores[, round_i] <- Matrix::colMeans(pseudo_norm[sampled_genes, , drop = FALSE])
  }

  background_score <- rowMeans(background_scores)
  cilia_score <- as.numeric(target_score) - background_score

  data.frame(
    pseudo_id = colnames(pseudo_norm),
    CiliaScore = cilia_score,
    target_score = as.numeric(target_score),
    background_score = as.numeric(background_score),
    stringsAsFactors = FALSE
  )
}

summarize_scores_by_sample <- function(score) {
  aggregate_one <- function(x) {
    c(
      mean_CiliaScore = mean(x$CiliaScore, na.rm = TRUE),
      median_CiliaScore = median(x$CiliaScore, na.rm = TRUE),
      sd_CiliaScore = stats::sd(x$CiliaScore, na.rm = TRUE),
      n_pseudocells = nrow(x),
      total_cells = sum(x$n_cells, na.rm = TRUE)
    )
  }

  split_cols <- c("sample_id", "group", "cell_type", "subcluster")
  key <- do.call(paste, c(score[split_cols], sep = "||"))
  chunks <- split(score, key)
  out <- do.call(rbind, lapply(chunks, aggregate_one))
  meta <- do.call(rbind, lapply(chunks, function(x) x[1, split_cols, drop = FALSE]))
  out <- cbind(meta, as.data.frame(out, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out$n_pseudocells <- as.integer(out$n_pseudocells)
  out$total_cells <- as.integer(out$total_cells)
  out
}

save_pseudoCiliaScore_results <- function(result, outdir) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  write.csv(result$score, file.path(outdir, "score.csv"), row.names = FALSE)
  write.csv(result$sample_summary, file.path(outdir, "sample_summary.csv"), row.names = FALSE)
  write.csv(result$gene_info, file.path(outdir, "gene_info.csv"), row.names = FALSE)
  writeLines(result$cilia_genes_used, file.path(outdir, "cilia_genes_used.txt"))
  writeLines(result$cilia_genes_missing, file.path(outdir, "cilia_genes_missing.txt"))
  saveRDS(result$background_mapping, file.path(outdir, "background_mapping.rds"))
  saveRDS(result, file.path(outdir, "result.rds"))
  invisible(outdir)
}

.match_cilia_genes <- function(cilia_genes,
                               matrix_genes,
                               case_insensitive_genes = FALSE) {
  cilia_genes <- unique(as.character(cilia_genes))
  cilia_genes <- cilia_genes[!is.na(cilia_genes) & cilia_genes != ""]

  if (!case_insensitive_genes) {
    used <- intersect(cilia_genes, matrix_genes)
    missing <- setdiff(cilia_genes, used)
    return(list(used = used, missing = missing))
  }

  gene_lookup <- setNames(matrix_genes, toupper(matrix_genes))
  requested <- toupper(cilia_genes)
  matched <- unname(gene_lookup[requested])
  used <- unique(matched[!is.na(matched)])
  missing <- cilia_genes[is.na(matched)]
  list(used = used, missing = missing)
}

cal_pseudoCiliaScore <- function(
  object,
  cilia_genes,
  meta = NULL,
  assay = "RNA",
  layer = "counts",
  sample_col = "sample_id",
  group_col = "group",
  celltype_col = NULL,
  subcluster_col = NULL,
  target_celltype = NULL,
  pseudo_mode = c("chunk", "sample"),
  cells_per_pseudo = 80,
  min_cells_per_pseudo = 30,
  min_total_count = 10,
  min_detect_pseudocells = 3,
  min_detect_rate = 0.01,
  remove_mt = TRUE,
  remove_ribo = FALSE,
  scale.factor = 10000,
  n.mean.bins = 20,
  n.detect.bins = 10,
  n.control = 50,
  num.rounds = 1000,
  seed = 123,
  case_insensitive_genes = FALSE,
  return_matrices = FALSE,
  verbose = TRUE
) {
  pseudo_mode <- match.arg(pseudo_mode)

  extracted <- extract_counts_and_meta(
    object = object,
    meta = meta,
    assay = assay,
    layer = layer
  )
  counts <- extracted$counts
  meta <- extracted$meta

  subsetted <- subset_target_cells(
    counts = counts,
    meta = meta,
    celltype_col = celltype_col,
    target_celltype = target_celltype
  )
  counts <- subsetted$counts
  meta <- subsetted$meta

  gene_match_initial <- .match_cilia_genes(
    cilia_genes = cilia_genes,
    matrix_genes = rownames(counts),
    case_insensitive_genes = case_insensitive_genes
  )
  if (length(gene_match_initial$used) < 5L) {
    stop("Too few cilia genes overlap with raw counts before filtering: ",
         length(gene_match_initial$used), ". At least 5 are required.")
  }

  pseudo <- build_pseudocells(
    counts = counts,
    meta = meta,
    sample_col = sample_col,
    group_col = group_col,
    celltype_col = celltype_col,
    subcluster_col = subcluster_col,
    pseudo_mode = pseudo_mode,
    cells_per_pseudo = cells_per_pseudo,
    min_cells_per_pseudo = min_cells_per_pseudo,
    seed = seed,
    verbose = verbose
  )

  pseudo_counts_filtered <- filter_genes_for_scoring(
    pseudo_counts = pseudo$pseudo_counts,
    cilia_genes = gene_match_initial$used,
    min_total_count = min_total_count,
    min_detect_pseudocells = min_detect_pseudocells,
    min_detect_rate = min_detect_rate,
    remove_mt = remove_mt,
    remove_ribo = remove_ribo
  )

  gene_match_filtered <- .match_cilia_genes(
    cilia_genes = cilia_genes,
    matrix_genes = rownames(pseudo_counts_filtered),
    case_insensitive_genes = case_insensitive_genes
  )
  cilia_genes_used <- gene_match_filtered$used
  cilia_genes_missing <- setdiff(unique(as.character(cilia_genes)), cilia_genes_used)

  if (length(cilia_genes_used) < 5L) {
    stop("Too few cilia genes remain after pseudo-cell gene filtering: ",
         length(cilia_genes_used), ". At least 5 are required.")
  }

  pseudo_norm <- normalize_pseudocounts(
    pseudo_counts = pseudo_counts_filtered,
    scale.factor = scale.factor
  )

  gene_info <- compute_gene_background_features(
    pseudo_counts = pseudo_counts_filtered,
    pseudo_norm = pseudo_norm,
    cilia_genes_used = cilia_genes_used,
    n.mean.bins = n.mean.bins,
    n.detect.bins = n.detect.bins
  )

  background_mapping <- match_background_genes(
    gene_info = gene_info,
    cilia_genes_used = cilia_genes_used,
    n.control = n.control,
    verbose = verbose
  )

  score_core <- calculate_score_with_random_controls(
    pseudo_norm = pseudo_norm,
    cilia_genes_used = cilia_genes_used,
    background_mapping = background_mapping,
    n.control = n.control,
    num.rounds = num.rounds,
    seed = seed,
    verbose = verbose
  )

  score <- merge(score_core, pseudo$pseudo_meta, by = "pseudo_id", all.x = TRUE, sort = FALSE)
  score <- score[match(score_core$pseudo_id, score$pseudo_id), , drop = FALSE]
  sample_summary <- summarize_scores_by_sample(score)

  parameters <- list(
    assay = assay,
    layer = layer,
    sample_col = sample_col,
    group_col = group_col,
    celltype_col = celltype_col,
    subcluster_col = subcluster_col,
    target_celltype = target_celltype,
    pseudo_mode = pseudo_mode,
    cells_per_pseudo = cells_per_pseudo,
    min_cells_per_pseudo = min_cells_per_pseudo,
    min_total_count = min_total_count,
    min_detect_pseudocells = min_detect_pseudocells,
    min_detect_rate = min_detect_rate,
    remove_mt = remove_mt,
    remove_ribo = remove_ribo,
    scale.factor = scale.factor,
    n.mean.bins = n.mean.bins,
    n.detect.bins = n.detect.bins,
    n.control = n.control,
    num.rounds = num.rounds,
    seed = seed,
    case_insensitive_genes = case_insensitive_genes,
    return_matrices = return_matrices
  )

  result <- list(
    score = score,
    sample_summary = sample_summary,
    pseudo_meta = pseudo$pseudo_meta,
    pseudo_counts = if (return_matrices) pseudo_counts_filtered else NULL,
    pseudo_norm = if (return_matrices) pseudo_norm else NULL,
    gene_info = gene_info,
    cilia_genes_used = cilia_genes_used,
    cilia_genes_missing = cilia_genes_missing,
    background_mapping = background_mapping,
    parameters = parameters
  )

  class(result) <- c("pseudoCiliaScore", class(result))
  result
}
