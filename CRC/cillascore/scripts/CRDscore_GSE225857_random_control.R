library(ggplot2)
library(ggpubr)
library(Seurat)
library(CRDscore)
library(qs)
library(Matrix)
library(readr)

prepare_crdscore_input <- function(expr, n.bins = 50) {
  expr <- as.matrix(expr)
  genes <- rownames(expr)
  genes.mean <- rowMeans(expr)
  expr.scaled <- sweep(expr, 1, genes.mean)
  dist <- 10 * ((2^expr) - 1)
  genes.dist <- log2(rowMeans(dist, na.rm = TRUE) + 1)
  genes.dist.factor <- arules::discretize(genes.dist, method = "frequency", breaks = n.bins)
  genes.dist.bins <- genes.dist.factor |>
    plyr::mapvalues(levels(genes.dist.factor), seq_along(levels(genes.dist.factor))) |>
    as.numeric() |>
    as.matrix()

  list(
    expr = expr,
    genes = genes,
    genes.mean = genes.mean,
    expr.scaled = expr.scaled,
    genes.dist = genes.dist,
    genes.dist.bins = genes.dist.bins
  )
}

calculate_crdscore_from_input <- function(crd_input, genes, seed = TRUE) {
  b.sign <- is.element(crd_input$genes, genes)
  if (sum(b.sign) <= 5) {
    stop("Non-enough-overlapping genes to calculate CRDscore")
  }
  if (seed) {
    set.seed(123)
  }
  r.scores <- CRDscore:::get_random_CRD(
    crd_input,
    crd_input$genes.dist.bins,
    b.sign,
    num.rounds = 1000
  )
  raw.scores <- colMeans(crd_input$expr.scaled[b.sign, ])

  # Keep the original CRDscore package direction: background score - target score.
  r.scores - raw.scores
}

# --- 配置 ---
WORK_DIR <- '/mnt/disk1/qiuzerui/downloads/CRC/cillascore/GSE225857'
DATASET <- basename(WORK_DIR)
FILES_DIR <- file.path(WORK_DIR, 'files')
PLOTS_DIR <- file.path(WORK_DIR, 'plots')
dir.create(FILES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

seed_base <- 12345
n_random_sets <- 100
random_run_label <- paste0("random_", n_random_sets)
set.seed(seed_base)

# --- 读取 gene set ---
gene_csv <- "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table1_primary_cilia_single_gene_confirmed.csv"
gene_set_name <- tools::file_path_sans_ext(basename(gene_csv))
gene_set_prefix <- gsub("[^A-Za-z0-9_]+", "_", gene_set_name)
gene_set_prefix <- gsub("_+", "_", gene_set_prefix)
gene_set_prefix <- gsub("^_|_$", "", gene_set_prefix)
gs <- read.csv(gene_csv, header = TRUE, check.names = FALSE)
gene_col <- intersect(c("gene", "Gene", "symbol", "Symbol", "gene_symbol", "GeneSymbol"), colnames(gs))[1]
if (is.na(gene_col)) {
  stop("No gene symbol column found in gene set CSV.")
}
target_genes_all <- as.character(gs[[gene_col]])
target_genes_all <- unique(target_genes_all[!is.na(target_genes_all) & target_genes_all != ""])
cat(sprintf("Gene set %s: %d genes\n", gene_set_name, length(target_genes_all)))

# --- 读取 Seurat 对象 ---
pbmc1 <- qread('/mnt/disk1/qiuzerui/downloads/CRC/GSE225857/Malignant.qs')
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 1000000)

# --- 过滤全 0 基因 ---
pbmc1 <- pbmc1[Matrix::rowSums(GetAssayData(pbmc1, layer = "counts")) > 0, ]

# --- 提取表达矩阵并转 log2 ---
seurat_data <- LayerData(pbmc1, assay = "RNA", layer = "data")
data_log2 <- as.matrix(seurat_data / log(2))
crd_input <- prepare_crdscore_input(data_log2, n.bins = 50)

# --- 匹配基因 ---
target_genes <- intersect(target_genes_all, rownames(pbmc1))
cat(sprintf("Matched %d / %d genes in expression matrix\n", length(target_genes), length(target_genes_all)))

# --- 合并 metadata 和临床信息 ---
meta <- pbmc1@meta.data
meta$id <- rownames(meta)
cli <- read.csv('/mnt/disk1/qiuzerui/downloads/CRC/GSE225857/GSE225857_Cli.csv', header = TRUE, check.names = FALSE)
rt <- merge(meta, cli, by.x = "orig.ident", by.y = "Patient ID")
rt$group <- ifelse(rt[['Liver metastasis (LM)']] == "Yes", "Metastasis", "Primary")

# --- Observed CRDscore ---
cat("Calculating observed CRDscore...\n")
score_obs <- calculate_crdscore_from_input(crd_input, target_genes, seed = TRUE)
gc()
score_obs <- as.numeric(score_obs)
names(score_obs) <- colnames(data_log2)

# --- 构建 score 表 ---
obs_df <- data.frame(cell_id = colnames(data_log2), CRDscore = score_obs, stringsAsFactors = FALSE)
obs_df <- merge(obs_df, rt[, c("id", "orig.ident", "group")], by.x = "cell_id", by.y = "id")
colnames(obs_df)[colnames(obs_df) == "orig.ident"] <- "sample_id"

write_csv(obs_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_observed_score.csv")))

# --- Observed 统计 ---
obs_primary <- obs_df$CRDscore[obs_df$group == "Primary"]
obs_meta <- obs_df$CRDscore[obs_df$group == "Metastasis"]
obs_delta <- mean(obs_meta) - mean(obs_primary)
obs_wt <- wilcox.test(CRDscore ~ group, data = obs_df)
obs_p <- obs_wt$p.value
obs_sd <- sqrt(((length(obs_primary)-1)*sd(obs_primary)^2 + (length(obs_meta)-1)*sd(obs_meta)^2) / (length(obs_primary)+length(obs_meta)-2))
obs_eff <- if (obs_sd > 0) obs_delta / obs_sd else NA
obs_rank_sum <- sum(rank(obs_df$CRDscore)[obs_df$group == "Metastasis"])
obs_auc <- (obs_rank_sum - length(obs_meta)*(length(obs_meta)+1)/2) / (length(obs_meta)*length(obs_primary))

cat(sprintf("Observed: delta=%.4f, p=%.2e, AUC=%.3f, effect_size=%.4f\n", obs_delta, obs_p, obs_auc, obs_eff))

# --- Random CRDscore null distribution ---
cat(sprintf("Calculating %d random CRDscore replicates...\n", n_random_sets))
random_candidate_genes <- setdiff(rownames(data_log2), target_genes_all)
if (length(random_candidate_genes) < length(target_genes)) {
  stop("Random candidate gene pool is smaller than the matched target gene set.")
}

cell_group <- rt$group[match(colnames(data_log2), rt$id)]
if (any(is.na(cell_group))) {
  stop("Some expression matrix cells do not have Primary/Metastasis group labels.")
}

score_sum <- numeric(length(colnames(data_log2)))
names(score_sum) <- colnames(data_log2)
random_gene_long <- vector("list", n_random_sets)
random_stats_list <- vector("list", n_random_sets)

for (random_id in seq_len(n_random_sets)) {
  set.seed(seed_base + random_id)
  rand_genes <- sample(random_candidate_genes, length(target_genes), replace = FALSE)

  score_rand <- calculate_crdscore_from_input(crd_input, rand_genes, seed = TRUE)
  score_rand <- as.numeric(score_rand)
  names(score_rand) <- colnames(data_log2)
  score_sum <- score_sum + score_rand

  rand_primary <- score_rand[cell_group == "Primary"]
  rand_meta <- score_rand[cell_group == "Metastasis"]
  rand_delta <- mean(rand_meta) - mean(rand_primary)
  rand_p <- wilcox.test(score_rand ~ cell_group)$p.value
  rand_sd <- sqrt(((length(rand_primary)-1)*sd(rand_primary)^2 + (length(rand_meta)-1)*sd(rand_meta)^2) / (length(rand_primary)+length(rand_meta)-2))
  rand_eff <- if (rand_sd > 0) rand_delta / rand_sd else NA
  rand_rank_sum <- sum(rank(score_rand)[cell_group == "Metastasis"])
  rand_auc <- (rand_rank_sum - length(rand_meta)*(length(rand_meta)+1)/2) / (length(rand_meta)*length(rand_primary))

  random_gene_long[[random_id]] <- data.frame(
    random_id = paste0("random_", random_id),
    gene = rand_genes,
    stringsAsFactors = FALSE
  )
  random_stats_list[[random_id]] <- data.frame(
    group = paste0("random_", random_id),
    n_genes = length(rand_genes),
    delta_mean = rand_delta,
    p_value = rand_p,
    neg_log10_p = -log10(pmax(rand_p, .Machine$double.xmin)),
    effect_size = rand_eff,
    AUC = rand_auc,
    direction = ifelse(rand_delta > 0, "Metastasis_high", "Primary_high"),
    stringsAsFactors = FALSE
  )

  if (random_id %% 50 == 0) {
    cat(sprintf("Finished %d / %d random signatures\n", random_id, n_random_sets))
    gc()
  }
}

random_stats_df <- do.call(rbind, random_stats_list)
random_genes_df <- do.call(rbind, random_gene_long)
write_csv(random_genes_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_", random_run_label, "_genes_long.csv")))

random_mean_df <- data.frame(
  cell_id = names(score_sum),
  CRDscore = as.numeric(score_sum / n_random_sets),
  stringsAsFactors = FALSE
)
random_mean_df <- merge(random_mean_df, rt[, c("id", "group")], by.x = "cell_id", by.y = "id")
write_csv(random_mean_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_mean_score.csv")))

# --- 保存统计结果 ---
observed_stats_df <- data.frame(
  group = "observed",
  n_genes = length(target_genes),
  delta_mean = obs_delta,
  p_value = obs_p,
  neg_log10_p = -log10(pmax(obs_p, .Machine$double.xmin)),
  effect_size = obs_eff,
  AUC = obs_auc,
  direction = ifelse(obs_delta > 0, "Metastasis_high", "Primary_high"),
  stringsAsFactors = FALSE
)
stats_df <- rbind(observed_stats_df, random_stats_df)
write_csv(stats_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_CRDscore_stats.csv")))
write_csv(random_stats_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_", random_run_label, "_stats.csv")))

# --- 绘图：Observed gene set CRDscore boxplot ---
obs_plot_df <- obs_df
obs_plot_df$group <- factor(obs_plot_df$group, levels = c("Primary", "Metastasis"))

p_observed <- ggplot(obs_plot_df, aes(x = group, y = CRDscore, fill = group)) +
  geom_violin(alpha = 0.35, width = 0.9, trim = TRUE) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.55) +
  theme_bw() + theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5)) +
  scale_fill_manual(values = c("Primary" = "#4DBBD5", "Metastasis" = "#E64B35")) +
  ylab("Gene set CRDscore") +
  ggtitle(sprintf("%s %s Gene set CRDscore\nDelta=%.4f, p=%.2e, AUC=%.3f",
                  DATASET, gene_set_name, obs_delta, obs_p, obs_auc))

ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_gene_set_CRDscore_boxplot.pdf")), plot = p_observed, width = 6, height = 5)
ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_gene_set_CRDscore_boxplot.png")), plot = p_observed, width = 6, height = 5, dpi = 300)

# --- 绘图：Random mean CRDscore boxplot ---
random_mean_plot_df <- random_mean_df
random_mean_plot_df$group <- factor(random_mean_plot_df$group, levels = c("Primary", "Metastasis"))
random_mean_primary <- random_mean_plot_df$CRDscore[random_mean_plot_df$group == "Primary"]
random_mean_meta <- random_mean_plot_df$CRDscore[random_mean_plot_df$group == "Metastasis"]
random_mean_delta <- mean(random_mean_meta) - mean(random_mean_primary)

p_random_mean <- ggplot(random_mean_plot_df, aes(x = group, y = CRDscore, fill = group)) +
  geom_violin(alpha = 0.35, width = 0.9, trim = TRUE) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.55) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 14, face = "bold", color = "black"),
        axis.text.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.y = element_text(size = 14, face = "bold", color = "black"),
        axis.title.x = element_blank(),
        legend.position = "none",
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5)) +
  scale_fill_manual(values = c("Primary" = "#4DBBD5", "Metastasis" = "#E64B35")) +
  ylab("Random mean CRDscore") +
  ggtitle(sprintf("%s %s Random mean CRDscore\nDelta=%.4f, n=%d",
                  DATASET, gene_set_name, random_mean_delta, n_random_sets))

ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_mean_CRDscore_boxplot.pdf")), plot = p_random_mean, width = 6, height = 5)
ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_mean_CRDscore_boxplot.png")), plot = p_random_mean, width = 6, height = 5, dpi = 300)

# --- 绘图：Random null distributions ---
plot_null_distribution <- function(random_df, metric, observed_value, xlab, output_prefix) {
  p <- ggplot(random_df, aes(x = .data[[metric]])) +
    geom_histogram(aes(y = after_stat(density)), bins = 40, fill = "grey75", color = "white") +
    geom_density(color = "grey30", linewidth = 0.8) +
    geom_vline(xintercept = observed_value, color = "#E64B35", linewidth = 1.1) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(size = 12, color = "black"),
          axis.title = element_text(size = 13, face = "bold", color = "black"),
          plot.title = element_text(face = "bold", size = 13, hjust = 0.5)) +
    xlab(xlab) +
    ylab("Density") +
    ggtitle(sprintf("%s %s random null distribution\nred line = observed", DATASET, gene_set_name))

  ggsave(paste0(output_prefix, ".pdf"), plot = p, width = 6, height = 5)
  ggsave(paste0(output_prefix, ".png"), plot = p, width = 6, height = 5, dpi = 300)
}

plot_null_distribution(
  random_stats_df,
  "delta_mean",
  obs_delta,
  "Random delta mean (Metastasis - Primary CRDscore)",
  file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_delta_mean_null_distribution"))
)
plot_null_distribution(
  random_stats_df,
  "AUC",
  obs_auc,
  "Random AUC",
  file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_AUC_null_distribution"))
)
plot_null_distribution(
  random_stats_df,
  "effect_size",
  obs_eff,
  "Random effect size",
  file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_effect_size_null_distribution"))
)
plot_null_distribution(
  random_stats_df,
  "neg_log10_p",
  -log10(pmax(obs_p, .Machine$double.xmin)),
  "Random -log10(p value)",
  file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_neg_log10_p_null_distribution"))
)

cat("\nDone! Results in:", WORK_DIR, "\n")
