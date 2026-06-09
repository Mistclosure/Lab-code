suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
})

source("/home/zerui/code/pseudo_cilia_score/R/pseudo_cilia_score.R")

PSEUDO_CILIA_GENE_FILES <- c(
  table2 = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table2_primary_cilia_single_gene_confirmed.csv",
  table1 = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table1_primary_cilia_single_gene_confirmed.csv"
)

read_pseudo_cilia_gene_list <- function(path) {
  df <- read.csv(path, header = TRUE, check.names = FALSE, fileEncoding = "UTF-8-BOM")
  gene_col <- grep("gene", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) {
    gene_col <- colnames(df)[1]
  }
  unique(as.character(df[[gene_col]]))
}

plot_pseudo_cilia_box <- function(data, x_col, y_col, title, ylab, out_png, sample_level = FALSE) {
  data <- data[!is.na(data[[x_col]]) & !is.na(data[[y_col]]), , drop = FALSE]
  if (nrow(data) == 0 || length(unique(data[[x_col]])) < 1) {
    warning("Skip empty plot: ", title)
    return(invisible(NULL))
  }
  data[[x_col]] <- factor(data[[x_col]], levels = unique(as.character(data[[x_col]])))
  groups <- levels(data[[x_col]])
  y_range <- range(data[[y_col]], na.rm = TRUE)
  y_span <- diff(y_range)
  if (!is.finite(y_span) || y_span == 0) {
    y_span <- max(abs(y_range), 1, na.rm = TRUE) * 0.1
  }
  y_label <- y_range[2] + y_span * 0.12
  y_limit <- y_range[2] + y_span * 0.28

  p <- ggplot(data, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[x_col]])) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.65) +
    geom_jitter(width = 0.15, size = if (sample_level) 2.2 else 0.9, alpha = if (sample_level) 0.9 else 0.35) +
    coord_cartesian(ylim = c(y_range[1], y_limit), clip = "off") +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(title = title, x = NULL, y = ylab)

  if (length(groups) == 2 && nrow(data) >= 2) {
    p <- p + stat_compare_means(
      comparisons = list(groups),
      method = "wilcox.test",
      label = "p.format",
      label.y = y_label,
      tip.length = 0.02
    )
  } else if (length(groups) > 2 && nrow(data) >= length(groups)) {
    p <- p + stat_compare_means(
      method = "kruskal.test",
      label = "p.format",
      label.y = y_label
    )
  }

  ggsave(out_png, p, width = 6, height = 5, dpi = 300)
  invisible(p)
}

write_pseudo_cilia_group_test <- function(sample_summary, group_col, out_file) {
  x <- sample_summary[!is.na(sample_summary[[group_col]]) & !is.na(sample_summary$mean_CiliaScore), , drop = FALSE]
  groups <- unique(x[[group_col]])
  if (length(groups) != 2 || nrow(x) < 4) {
    writeLines("Group test skipped: need exactly two groups and at least four sample-level rows.", out_file)
    return(invisible(NULL))
  }
  wt <- wilcox.test(mean_CiliaScore ~ .data_group, data = data.frame(
    mean_CiliaScore = x$mean_CiliaScore,
    .data_group = x[[group_col]]
  ))
  writeLines(c(
    paste0("Wilcoxon rank-sum test by ", group_col),
    paste0("groups = ", paste(groups, collapse = ", ")),
    paste0("p.value = ", signif(wt$p.value, 5)),
    "",
    "Note: pseudo-cell chunks are visualized but formal comparison should use sample-level summaries."
  ), out_file)
}

save_standard_pseudo_cilia_plots <- function(result,
                                            dataset,
                                            signature_name,
                                            outdir,
                                            group_col = "group",
                                            stage_col = NULL) {
  plot_dir <- file.path(outdir, "plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  plot_pseudo_cilia_box(
    result$score, group_col, "CiliaScore",
    paste0(signature_name, "+", dataset, "+pseudoCiliaScore (pseudo-cell ", group_col, ")"),
    "Pseudo-cell CiliaScore",
    file.path(plot_dir, paste0(signature_name, "_pseudo_cell_", group_col, ".png"))
  )

  plot_pseudo_cilia_box(
    result$sample_summary, group_col, "mean_CiliaScore",
    paste0(signature_name, "+", dataset, "+pseudoCiliaScore (sample ", group_col, ")"),
    "Sample mean CiliaScore",
    file.path(plot_dir, paste0(signature_name, "_sample_", group_col, ".png")),
    sample_level = TRUE
  )

  if (!is.null(stage_col) && stage_col %in% colnames(result$sample_summary)) {
    plot_pseudo_cilia_box(
      result$sample_summary, stage_col, "mean_CiliaScore",
      paste0(signature_name, "+", dataset, "+pseudoCiliaScore (sample ", stage_col, ")"),
      "Sample mean CiliaScore",
      file.path(plot_dir, paste0(signature_name, "_sample_", stage_col, ".png")),
      sample_level = TRUE
    )
  }

  p_scatter <- ggplot(result$score, aes(x = background_score, y = target_score, color = .data[[group_col]])) +
    geom_point(size = 1.2, alpha = 0.65) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0(signature_name, "+", dataset, "+target vs background"),
         x = "Background score", y = "Target score")
  ggsave(file.path(plot_dir, paste0(signature_name, "_target_vs_background.png")), p_scatter, width = 6, height = 5, dpi = 300)

  used_gene_info <- result$gene_info[result$gene_info$used_as_cilia_gene, , drop = FALSE]
  p_detect <- ggplot(used_gene_info, aes(x = detection_rate)) +
    geom_histogram(bins = 20, fill = "#4C78A8", color = "white") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0(signature_name, "+", dataset, "+used cilia gene detection"),
         x = "Detection rate in pseudo-cells", y = "Gene count")
  ggsave(file.path(plot_dir, paste0(signature_name, "_used_gene_detection_rate.png")), p_detect, width = 6, height = 5, dpi = 300)
}

run_pseudo_cilia_for_gene_sets <- function(obj,
                                           dataset,
                                           out_root,
                                           sample_col,
                                           group_col,
                                           stage_df = NULL,
                                           stage_colname = "Stage",
                                           cells_per_pseudo = 80,
                                           min_cells_per_pseudo = 30) {
  for (set_name in names(PSEUDO_CILIA_GENE_FILES)) {
    signature_file <- PSEUDO_CILIA_GENE_FILES[[set_name]]
    signature_name <- tools::file_path_sans_ext(basename(signature_file))
    cilia_genes <- read_pseudo_cilia_gene_list(signature_file)
    outdir <- file.path(out_root, signature_name)
    dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

    message("Running ", dataset, " pseudoCiliaScore: ", signature_name)
    result <- cal_pseudoCiliaScore(
      object = obj,
      cilia_genes = cilia_genes,
      assay = "RNA",
      layer = "counts",
      sample_col = sample_col,
      group_col = group_col,
      celltype_col = NULL,
      subcluster_col = NULL,
      target_celltype = NULL,
      pseudo_mode = "chunk",
      cells_per_pseudo = cells_per_pseudo,
      min_cells_per_pseudo = min_cells_per_pseudo,
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
    )

    if (!is.null(stage_df)) {
      result$score <- merge(result$score, stage_df, by = "sample_id", all.x = TRUE, sort = FALSE)
      result$sample_summary <- merge(result$sample_summary, stage_df, by = "sample_id", all.x = TRUE, sort = FALSE)
    }

    save_pseudoCiliaScore_results(result, outdir)
    save_standard_pseudo_cilia_plots(
      result = result,
      dataset = dataset,
      signature_name = signature_name,
      outdir = outdir,
      group_col = "group",
      stage_col = if (!is.null(stage_df)) stage_colname else NULL
    )
    write_pseudo_cilia_group_test(result$sample_summary, "group", file.path(outdir, "sample_group_wilcox.txt"))
  }
}
