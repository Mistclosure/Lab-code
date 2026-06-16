library(ggplot2)
library(ggpubr)
library(Seurat)
library(CRDscore)
library(qs)
library(Matrix)
library(readr)

# --- 配置 ---
WORK_DIR <- '/mnt/disk1/qiuzerui/downloads/CRC/cillascore/GSE132465'
DATASET <- basename(WORK_DIR)
FILES_DIR <- file.path(WORK_DIR, 'files')
PLOTS_DIR <- file.path(WORK_DIR, 'plots')
dir.create(FILES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOTS_DIR, recursive = TRUE, showWarnings = FALSE)

seed_base <- 12345
n_random_sets <- 3
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
pbmc1 <- qread('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/qs/Seurat/Malignant_RNA_assay.qs')
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 1000000)

# --- 提取表达矩阵并转 log2 ---
seurat_data <- LayerData(pbmc1, assay = "RNA", layer = "data")
data_log2 <- as.data.frame(seurat_data / log(2))

# --- 匹配基因 ---
target_genes <- intersect(target_genes_all, rownames(pbmc1))
cat(sprintf("Matched %d / %d genes in expression matrix\n", length(target_genes), length(target_genes_all)))

# --- 合并 metadata 和临床信息 ---
meta <- pbmc1@meta.data
meta$id <- rownames(meta)
cli <- read.csv('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/metadata/GSE132465_Cli.csv', header = TRUE, check.names = FALSE)
rt <- merge(meta, cli, by.x = "orig.ident", by.y = "Tumor")
rt$group <- ifelse(grepl("M0$", rt[['TNM stage']]), "Primary", "Metastasis")

# --- Observed CRDscore ---
cat("Calculating observed CRDscore...\n")
score_obs <- cal_CRDscore(expr = data_log2, n.bins = 50, circadians = target_genes, study.type = "scRNAseq")
gc()
score_obs <- as.numeric(score_obs)
names(score_obs) <- colnames(data_log2)

# --- 构建 score 表 ---
obs_df <- data.frame(cell_id = colnames(data_log2), CRDscore = score_obs, stringsAsFactors = FALSE)
obs_df <- merge(obs_df, rt[, c("id", "orig.ident", "group")], by.x = "cell_id", by.y = "id")
colnames(obs_df)[colnames(obs_df) == "orig.ident"] <- "sample_id"

write_csv(obs_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_observed_score.csv")))

# --- Random CRDscore ---
cat(sprintf("Calculating %d random CRDscore replicates...\n", n_random_sets))
random_candidate_genes <- setdiff(rownames(data_log2), target_genes_all)
if (length(random_candidate_genes) < length(target_genes)) {
  stop("Random candidate gene pool is smaller than the matched target gene set.")
}
random_score_dfs <- vector("list", n_random_sets)
random_gene_sets <- vector("list", n_random_sets)

for (random_id in seq_len(n_random_sets)) {
  set.seed(seed_base + random_id)
  rand_genes <- sample(random_candidate_genes, length(target_genes), replace = FALSE)
  random_gene_sets[[random_id]] <- rand_genes

  score_rand <- cal_CRDscore(expr = data_log2, n.bins = 50, circadians = rand_genes, study.type = "scRNAseq")
  gc()
  score_rand <- as.numeric(score_rand)
  names(score_rand) <- colnames(data_log2)

  rand_df <- data.frame(cell_id = colnames(data_log2), CRDscore = score_rand, stringsAsFactors = FALSE)
  rand_df <- merge(rand_df, rt[, c("id", "group")], by.x = "cell_id", by.y = "id")
  rand_df$random_id <- paste0("random_", random_id)
  random_score_dfs[[random_id]] <- rand_df

  write_csv(
    data.frame(random_id = rand_df$random_id[1], gene = rand_genes, stringsAsFactors = FALSE),
    file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_", random_id, "_genes.csv"))
  )
  write_csv(rand_df, file.path(FILES_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_", random_id, "_score.csv")))
}

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

random_stats_df <- do.call(rbind, lapply(seq_len(n_random_sets), function(random_id) {
  rand_df <- random_score_dfs[[random_id]]
  rand_primary <- rand_df$CRDscore[rand_df$group == "Primary"]
  rand_meta <- rand_df$CRDscore[rand_df$group == "Metastasis"]
  rand_delta <- mean(rand_meta) - mean(rand_primary)
  rand_wt <- wilcox.test(CRDscore ~ group, data = rand_df)
  rand_p <- rand_wt$p.value
  rand_sd <- sqrt(((length(rand_primary)-1)*sd(rand_primary)^2 + (length(rand_meta)-1)*sd(rand_meta)^2) / (length(rand_primary)+length(rand_meta)-2))
  rand_eff <- if (rand_sd > 0) rand_delta / rand_sd else NA
  rand_rank_sum <- sum(rank(rand_df$CRDscore)[rand_df$group == "Metastasis"])
  rand_auc <- (rand_rank_sum - length(rand_meta)*(length(rand_meta)+1)/2) / (length(rand_meta)*length(rand_primary))

  cat(sprintf("Random %d: delta=%.4f, p=%.2e, AUC=%.3f, effect_size=%.4f\n", random_id, rand_delta, rand_p, rand_auc, rand_eff))

  data.frame(
    group = paste0("random_", random_id),
    n_genes = length(random_gene_sets[[random_id]]),
    delta_mean = rand_delta,
    p_value = rand_p,
    neg_log10_p = -log10(pmax(rand_p, .Machine$double.xmin)),
    effect_size = rand_eff,
    AUC = rand_auc,
    direction = ifelse(rand_delta > 0, "Metastasis_high", "Primary_high"),
    stringsAsFactors = FALSE
  )
}))

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

# --- 绘图：Random CRDscore boxplots ---
for (random_id in seq_len(n_random_sets)) {
  rand_plot_df <- random_score_dfs[[random_id]]
  rand_plot_df$group <- factor(rand_plot_df$group, levels = c("Primary", "Metastasis"))
  rand_stats <- random_stats_df[random_stats_df$group == paste0("random_", random_id), ]

  p_random <- ggplot(rand_plot_df, aes(x = group, y = CRDscore, fill = group)) +
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
    ylab("Random CRDscore") +
    ggtitle(sprintf("%s %s Random %d CRDscore\nDelta=%.4f, p=%.2e, AUC=%.3f",
                    DATASET, gene_set_name, random_id, rand_stats$delta_mean, rand_stats$p_value, rand_stats$AUC))

  ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_", random_id, "_CRDscore_boxplot.pdf")), plot = p_random, width = 6, height = 5)
  ggsave(file.path(PLOTS_DIR, paste0(DATASET, "_", gene_set_prefix, "_random_", random_id, "_CRDscore_boxplot.png")), plot = p_random, width = 6, height = 5, dpi = 300)
}

cat("\nDone! Results in:", WORK_DIR, "\n")
