# ==============================================================================
# 脚本二：癌细胞功能分析与相关性计算 (02_Malignant_Analysis.R)
# 依赖：必须先运行 01_Tcell_Analysis.R，确保工作目录下有 
#       scRNA_annotated.qs 和 Sample_Level_Tcell_Proportion.csv
# ==============================================================================

# 加载必要的库
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(CRDscore)
library(qs)
library(Matrix)

setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE178341/')

# 全局分析参数 (需与脚本一保持一致)
sample_col <- "orig.ident"    
stage_col  <- "Tumor.Stage"   

# ------------------------------------------------------------------------------
# 1. 加载数据与准备
# ------------------------------------------------------------------------------
cat(">>> 正在加载基础对象与 T 细胞占比数据...\n")
scRNA <- qread('scRNA.qs')
sample_stats <- read.csv("Sample_Level_Tcell_Proportion.csv")

# 读取目标基因集
CRC_data <- read.csv("ciliopathy_genes.csv", header = T, check.names = F)
target_genes <- as.character(CRC_data[,1])
target_genes <- intersect(target_genes, rownames(scRNA)) 

# 提取癌细胞并获取表达矩阵
scRNA_malig <- subset(scRNA, is_malignant == TRUE)
qsave(scRNA_malig,'scRNA_malig.qs')
exp_malig <- LayerData(scRNA_malig, assay = "RNA", layer = "data_log2") 

# ------------------------------------------------------------------------------
# 2. 癌细胞 172 基因总体 CRDscore 与 T 细胞占比相关性
# ------------------------------------------------------------------------------
cat(">>> 正在为癌细胞计算 172 基因总体 CRDscore...\n")
malig_overall_scores <- cal_CRDscore(expr = as.data.frame(exp_malig), n.bins = 50, 
                                     circadians = target_genes, study.type = "scRNAseq")

scRNA_malig$malig_score <- NA
scRNA_malig@meta.data[names(malig_overall_scores), "malig_score"] <- malig_overall_scores

plot_data_malig_score <- scRNA_malig@meta.data %>%
  group_by(!!sym(sample_col), !!sym(stage_col)) %>%
  summarise(Avg_Malig_CRDScore = mean(malig_score, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(sample_stats %>% select(!!sym(sample_col), Tcell_Pct), by = sample_col)

p3 <- ggplot(plot_data_malig_score, aes(x = Avg_Malig_CRDScore, y = Tcell_Pct)) +
  geom_point(aes(color = !!sym(stage_col)), size = 4) +
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  stat_cor(method = "pearson") + 
  theme_bw() +
  labs(title = "Correlation: Malignant 172-gene Score vs. T Cell %",
       x = "Average Malignant Cell CRDscore (per Sample)",
       y = "T Cell Percentage per Sample (%)")

print(p3)
ggsave("p3_Malignant_Score_vs_T_cell_Percentage.png", p3, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3. 极速版：全基因组 平均表达量 vs T细胞占比 相关性计算
# ------------------------------------------------------------------------------
cat(">>> 正在启动矩阵优化版计算 (平均表达量) ...\n")
sample_ids <- scRNA_malig$orig.ident
mm <- model.matrix(~ 0 + sample_ids)
colnames(mm) <- levels(as.factor(sample_ids))

sample_mean_matrix <- (exp_malig %*% mm) %*% diag(1/colSums(mm))
colnames(sample_mean_matrix) <- colnames(mm)
sample_mean_matrix <- as.matrix(sample_mean_matrix)

common_samples <- intersect(colnames(sample_mean_matrix), sample_stats[[sample_col]])
t_pct_vector <- sample_stats$Tcell_Pct[match(common_samples, sample_stats[[sample_col]])]
expr_matrix_sub <- sample_mean_matrix[, common_samples]

cors <- apply(expr_matrix_sub, 1, function(x) {
  if(sd(x) == 0) return(c(NA, NA))
  res <- cor.test(x, t_pct_vector, method = "pearson")
  return(c(res$estimate, res$p.value))
})

all_gene_cor_results <- as.data.frame(t(cors))
colnames(all_gene_cor_results) <- c("Pearson_R", "P_Value")
all_gene_cor_results$Gene <- rownames(all_gene_cor_results)
all_gene_cor_results <- all_gene_cor_results %>%
  filter(!is.na(Pearson_R)) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  arrange(desc(Pearson_R))

write.csv(all_gene_cor_results, "Malignant_AllGenes_Expr_vs_Tcell_Pct_Fast.csv", row.names = FALSE)
target_gene_cor_results <- all_gene_cor_results %>% filter(Gene %in% target_genes)
write.csv(target_gene_cor_results, "Malignant_TargetGenes_Expr_vs_Tcell_Pct1.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. 极速版：全基因组 非零值比例 vs T细胞占比 相关性计算 (Step 6)
# ------------------------------------------------------------------------------
cat(">>> 正在启动极速版计算 (非零检出比例) ...\n")

exp_malig_bin <- exp_malig
exp_malig_bin@x <- rep(1, length(exp_malig_bin@x))

nonzero_counts <- exp_malig_bin %*% mm
nonzero_props <- nonzero_counts %*% diag(1/colSums(mm))
colnames(nonzero_props) <- colnames(mm)
nonzero_props <- as.matrix(nonzero_props)

prop_matrix_sub <- nonzero_props[, common_samples]

cors_prop <- apply(prop_matrix_sub, 1, function(x) {
  if(sd(x) == 0) return(c(NA, NA)) 
  res <- cor.test(x, t_pct_vector, method = "pearson")
  return(c(res$estimate, res$p.value))
})

prop_cor_results <- as.data.frame(t(cors_prop))
colnames(prop_cor_results) <- c("Pearson_R", "P_Value")
prop_cor_results$Gene <- rownames(prop_cor_results)

prop_cor_results <- prop_cor_results %>%
  filter(!is.na(Pearson_R)) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  arrange(desc(Pearson_R))

write.csv(prop_cor_results, "Malignant_AllGenes_NonzeroProp_vs_Tcell_Pct_Fast.csv", row.names = FALSE)

target_prop_cor_results <- prop_cor_results %>% filter(Gene %in% target_genes)
write.csv(target_prop_cor_results, "Malignant_TargetGenes_NonzeroProp_vs_Tcell_Pct.csv", row.names = FALSE)

cat(">>> 癌细胞相关性极速版分析全部完成！\n")