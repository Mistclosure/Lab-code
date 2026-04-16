# 加载必要的库
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(CRDscore)
library(SingleR) 
library(celldex) 
library(qs)

setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE178341/')

# ==============================================================================
# Part I: 数据加载、注释与预处理
# ==============================================================================

# 1. 加载基础数据与参考集
load('scRNA.Rdata')
ref <- HumanPrimaryCellAtlasData()

# 2. 数据标准化与 log2 转换 (为 SingleR 和 CRDscore 做准备)
scRNA <- NormalizeData(scRNA, assay = "RNA") 
scRNA[["RNA"]]$data_log2 <- GetAssayData(scRNA, assay = "RNA", layer = "data") / log(2)

# 3. 运行 SingleR 细胞自动注释
test_counts <- GetAssayData(scRNA, assay = "RNA", layer = "data")
pred <- SingleR(test = test_counts, ref = ref, labels = ref$label.main)
scRNA$SingleR <- pred$labels

# 4. 提取与合并临床元数据 (Cli.csv)
metadata <- scRNA@meta.data
metadata$cell_id_tmp <- rownames(metadata)
cli_data <- read.csv('Cli.csv') %>% distinct(PatientBarcode, .keep_all = TRUE)
metadata_merged <- left_join(metadata, cli_data, by = c("orig.ident" = "PatientBarcode"))
rownames(metadata_merged) <- metadata_merged$cell_id_tmp
metadata_merged$cell_id_tmp <- NULL 
scRNA@meta.data <- metadata_merged
# 所有底层数据处理完毕，保存对象
qsave(scRNA, 'scRNA.qs')
scRNA = qread('scRNA.qs')
# 5. 读取目标基因集
CRC_data <- read.csv("ciliopathy_genes.csv", header = T, check.names = F)
target_genes <- as.character(CRC_data[,1])
target_genes <- intersect(target_genes, rownames(scRNA)) 

# 6. 计算 T 细胞的总体 CRDscore   不计算step2，step3请跳过
scRNA_t <- subset(scRNA, cluster1 == "T_cell")
exp_t <- as.data.frame(LayerData(scRNA_t, assay = "RNA", layer = "data_log2"))

t_scores <- cal_CRDscore(expr = exp_t, n.bins = 50, circadians = target_genes, study.type = "scRNAseq")
scRNA$score <- NA
scRNA@meta.data[names(t_scores), "score"] <- t_scores

# 7. 标记癌细胞 (Malignant) 并保存最终对象
malig_df <- read.table('Malignant cells.txt', header = TRUE, stringsAsFactors = FALSE)
malignant_barcodes <- as.character(malig_df[, 1])
scRNA$is_malignant <- colnames(scRNA) %in% malignant_barcodes


# ==============================================================================
# Part II: 下游统计分析与绘图
# ==============================================================================

# 设置全局分析参数
tcell_label_name <- "T_cell" 
sample_col <- "orig.ident"    
stage_col  <- "Tumor.Stage"   
score_col  <- "score"         
meta_data <- scRNA@meta.data

# ------------------------------------------------------------------------------
# 1. 计算 T 细胞总占比
# ------------------------------------------------------------------------------
total_cells <- ncol(scRNA)
tcell_count <- sum(scRNA$cluster1 == tcell_label_name, na.rm = TRUE)
total_pct   <- (tcell_count / total_cells) * 100

cat("--------------------------------------------------\n")
cat(paste0(">>> 总细胞数: ", total_cells, "\n"))
cat(paste0(">>> T 细胞总数: ", tcell_count, "\n"))
cat(paste0(">>> T 细胞总占比: ", round(total_pct, 2), "%\n"))
cat("--------------------------------------------------\n")

# ------------------------------------------------------------------------------
# 2. 计算并绘制：基于 Tumor Stage 分组的 T 细胞占比
# ------------------------------------------------------------------------------
stage_stats <- meta_data %>%
  filter(!is.na(!!sym(stage_col))) %>%
  mutate(!!sym(stage_col) := as.factor(!!sym(stage_col))) %>% 
  group_by(!!sym(stage_col)) %>%
  summarise(
    Total_in_Stage = n(),
    Tcell_in_Stage = sum(cluster1 == tcell_label_name, na.rm = TRUE),
    Percentage = (Tcell_in_Stage / Total_in_Stage) * 100
  )

p1 <- ggplot(stage_stats, aes(x = !!sym(stage_col), y = Percentage, fill = !!sym(stage_col))) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), vjust = -0.5, size = 5) +
  scale_fill_viridis_d(option = "D") +
  theme_minimal() +
  labs(title = "T Cell Proportion across Tumor Stages",
       x = "Tumor Stage", y = "T Cell Percentage (%)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("p1_T_cell_Proportion_by_Stage.png", p1, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 3. 计算并绘制：样本水平 T 细胞占比与 T 细胞 Score 的相关性
# ------------------------------------------------------------------------------
sample_stats <- meta_data %>%
  filter(!is.na(!!sym(stage_col))) %>%
  group_by(!!sym(sample_col), !!sym(stage_col)) %>%
  summarise(
    Tcell_Pct = (sum(cluster1 == tcell_label_name, na.rm = TRUE) / n()) * 100,
    Avg_Score = mean(!!sym(score_col), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(Avg_Score))

p2 <- ggplot(sample_stats, aes(x = Avg_Score, y = Tcell_Pct)) +
  geom_point(aes(color = !!sym(stage_col)), size = 4) +
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
  theme_bw() +
  labs(title = paste0("Correlation: T Cell % vs. T-Cell ", score_col),
       x = "Average T-Cell Score (per Sample)",
       y = "T Cell Percentage per Sample (%)")

print(p2)
ggsave("p2_T_cell_Score_vs_Percentage.png", p2, width = 8, height = 6, dpi = 300)

# ------------------------------------------------------------------------------
# 4. 计算并绘制：癌细胞 172 基因总体 CRDscore 与 T 细胞占比相关性
# ------------------------------------------------------------------------------
cat(">>> 正在为癌细胞计算 172 基因总体 CRDscore...\n")
# --- 极速版前置计算：计算 T 细胞占比并生成表格 ---
sample_stats <- scRNA@meta.data %>%
  filter(!is.na(!!sym(stage_col))) %>%
  group_by(!!sym(sample_col), !!sym(stage_col)) %>%
  summarise(
    Total_Cells = n(),                                     # 该样本总细胞数
    Tcell_Count = sum(cluster1 == tcell_label_name, na.rm = TRUE), # 该样本 T 细胞数
    Tcell_Pct = (Tcell_Count / Total_Cells) * 100,         # T 细胞占比
    .groups = 'drop'
  )

# 【新增】将 T 细胞占比输出为 CSV 表格
write.csv(sample_stats, "Sample_Level_Tcell_Proportion.csv", row.names = FALSE)
cat(">>> T 细胞占比结果已保存至: Sample_Level_Tcell_Proportion.csv\n")
scRNA_malig <- subset(scRNA, is_malignant == TRUE)
qsave(scRNA_malig,'scRNA_malig.qs')

exp_malig <- LayerData(scRNA_malig, assay = "RNA", layer = "data_log2") 

malig_overall_scores <- cal_CRDscore(expr = as.data.frame(exp_malig), n.bins = 50, 
                                     circadians = target_genes, study.type = "scRNAseq")

scRNA$malig_score <- NA
scRNA@meta.data[names(malig_overall_scores), "malig_score"] <- malig_overall_scores

plot_data_malig_score <- scRNA@meta.data %>%
  filter(is_malignant == TRUE) %>%
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


# ==============================================================================
# 5. 极速版：全基因组相关性计算 (优化矩阵运算)
# ==============================================================================
cat(">>> 正在启动矩阵优化版计算 (预计 2-5 分钟)...\n")

# 1. 预先计算所有基因在每个样本中的平均表达量 (这一步是核心优化)
# 使用 Seurat 的 AverageExpression 或底层矩阵乘法
sample_ids <- scRNA_malig$orig.ident
# 创建一个“细胞-样本”的对角矩阵（用于快速求均值）
mm <- model.matrix(~ 0 + sample_ids)
colnames(mm) <- levels(as.factor(sample_ids))

# 矩阵相乘直接获得所有基因在所有样本的表达总和，然后除以样本细胞数得到均值
# exp_malig 是 Genes x Cells, mm 是 Cells x Samples -> 结果是 Genes x Samples
sample_mean_matrix <- (exp_malig %*% mm) %*% diag(1/colSums(mm))
colnames(sample_mean_matrix) <- colnames(mm)
sample_mean_matrix <- as.matrix(sample_mean_matrix)

# 2. 准备 T 细胞占比向量，并确保样本顺序一致
common_samples <- intersect(colnames(sample_mean_matrix), sample_stats[[sample_col]])
t_pct_vector <- sample_stats$Tcell_Pct[match(common_samples, sample_stats[[sample_col]])]
expr_matrix_sub <- sample_mean_matrix[, common_samples]

# 3. 批量计算相关性 (使用最基础的 cor 函数，速度极快)
# cor() 可以直接计算矩阵每一行与一个向量的相关性
cors <- apply(expr_matrix_sub, 1, function(x) {
  if(sd(x) == 0) return(c(NA, NA))
  res <- cor.test(x, t_pct_vector, method = "pearson")
  return(c(res$estimate, res$p.value))
})

# 4. 整理结果
all_gene_cor_results <- as.data.frame(t(cors))
colnames(all_gene_cor_results) <- c("Pearson_R", "P_Value")
all_gene_cor_results$Gene <- rownames(all_gene_cor_results)

all_gene_cor_results <- all_gene_cor_results %>%
  filter(!is.na(Pearson_R)) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  arrange(desc(Pearson_R))

write.csv(all_gene_cor_results, "Malignant_AllGenes_Expr_vs_Tcell_Pct_Fast.csv", row.names = FALSE)
cat(">>> 极速版分析完成！\n")
target_gene_cor_results <- all_gene_cor_results %>%
  filter(Gene %in% target_genes)

write.csv(target_gene_cor_results, "Malignant_TargetGenes_Expr_vs_Tcell_Pct1.csv")
# ==============================================================================
# Step 6: 极速版 - 癌细胞中每个基因非零值比例与 T 细胞比例的相关性
# ==============================================================================
cat(">>> 正在启动 Step 6: 基因非零表达比例与 T细胞占比 相关性计算...\n")

# 1. 极速构建二值化稀疏矩阵 (巧妙利用 dgCMatrix 结构)
# 直接将稀疏矩阵中实际存储的非零值 (即 @x) 全部替换为 1
# 这一步是 O(N) 复杂度，比任何 >0 的逻辑判断运算都要快得多，且完全保留稀疏性
exp_malig_bin <- exp_malig
exp_malig_bin@x <- rep(1, length(exp_malig_bin@x))

# 2. 利用你在 Step 5 构建的 model.matrix (mm: 细胞 x 样本) 
# 通过矩阵乘法快速汇总：每个样本中，各个基因表达 > 0 的细胞总数
# 结果矩阵维度: Genes x Samples
nonzero_counts <- exp_malig_bin %*% mm

# 3. 计算非零细胞比例 (Proportion)
# colSums(mm) 刚好是每个样本中总癌细胞数，除一下就得到了比例
nonzero_props <- nonzero_counts %*% diag(1/colSums(mm))
colnames(nonzero_props) <- colnames(mm)
nonzero_props <- as.matrix(nonzero_props)

# 4. 对齐 T 细胞占比数据
common_samples_step7 <- intersect(colnames(nonzero_props), sample_stats[[sample_col]])
t_pct_vector_step7 <- sample_stats$Tcell_Pct[match(common_samples_step7, sample_stats[[sample_col]])]
prop_matrix_sub <- nonzero_props[, common_samples_step7]

# 5. 批量计算相关性 (非零比例 vs T细胞占比)
cat(">>> 正在计算非零比例与 T 细胞占比的相关性 (约耗时 1-3 分钟)...\n")
cors_prop <- apply(prop_matrix_sub, 1, function(x) {
  # 如果某个基因在所有样本里的非零比例都一样（方差为0），则无法计算相关性
  if(sd(x) == 0) return(c(NA, NA)) 
  res <- cor.test(x, t_pct_vector_step7, method = "pearson")
  return(c(res$estimate, res$p.value))
})

# 6. 整理并格式化结果表
prop_cor_results <- as.data.frame(t(cors_prop))
colnames(prop_cor_results) <- c("Pearson_R", "P_Value")
prop_cor_results$Gene <- rownames(prop_cor_results)

prop_cor_results <- prop_cor_results %>%
  filter(!is.na(Pearson_R)) %>%
  mutate(Adj_P_Value = p.adjust(P_Value, method = "fdr")) %>%
  arrange(desc(Pearson_R))

# 7. 导出全基因组与目标基因集的计算结果
write.csv(prop_cor_results, "Malignant_AllGenes_NonzeroProp_vs_Tcell_Pct_Fast.csv", row.names = FALSE)

target_prop_cor_results <- prop_cor_results %>%
  filter(Gene %in% target_genes)
write.csv(target_prop_cor_results, "Malignant_TargetGenes_NonzeroProp_vs_Tcell_Pct.csv", row.names = FALSE)

cat(">>> Step 6 极速版分析完成！非零比例相关性结果已成功保存。\n")