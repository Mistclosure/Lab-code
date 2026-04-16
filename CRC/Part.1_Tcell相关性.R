# ==============================================================================
# 脚本一：基础预处理与 T 细胞分析 (01_Tcell_Analysis.R)
# 请先运行此脚本，它会生成脚本二所需的 T 细胞占比数据和基础大对象
# ==============================================================================

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

# 全局分析参数
tcell_label_name <- "T_cell" 
sample_col <- "orig.ident"    
stage_col  <- "Tumor.Stage"   
score_col  <- "score"         

# ------------------------------------------------------------------------------
# Part I: 数据加载、注释与预处理
# ------------------------------------------------------------------------------
# 1. 加载基础数据与参考集
load('scRNA.Rdata')
ref <- HumanPrimaryCellAtlasData()

# 2. 数据标准化与 log2 转换
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

# 5. 读取目标基因集
CRC_data <- read.csv("ciliopathy_genes.csv", header = T, check.names = F)
target_genes <- as.character(CRC_data[,1])
target_genes <- intersect(target_genes, rownames(scRNA)) 

# 6. 标记癌细胞 (Malignant) - 提前标记，存入全局对象
malig_df <- read.table('Malignant cells.txt', header = TRUE, stringsAsFactors = FALSE)
malignant_barcodes <- as.character(malig_df[, 1])
scRNA$is_malignant <- colnames(scRNA) %in% malignant_barcodes
# 保存包含所有注释和 T 细胞打分的完整对象 (供脚本二读取)
qsave(scRNA, 'scRNA.qs')
# 7. 计算 T 细胞的总体 CRDscore 
scRNA_t <- subset(scRNA, cluster1 == "T_cell")
exp_t <- as.data.frame(LayerData(scRNA_t, assay = "RNA", layer = "data_log2"))
t_scores <- cal_CRDscore(expr = exp_t, n.bins = 50, circadians = target_genes, study.type = "scRNAseq")
scRNA$score <- NA
scRNA@meta.data[names(t_scores), "score"] <- t_scores


meta_data <- scRNA@meta.data

# ------------------------------------------------------------------------------
# Part II: T 细胞下游统计分析与绘图
# ------------------------------------------------------------------------------

# 1. 计算 T 细胞总占比
total_cells <- ncol(scRNA)
tcell_count <- sum(scRNA$cluster1 == tcell_label_name, na.rm = TRUE)
total_pct   <- (tcell_count / total_cells) * 100

cat("--------------------------------------------------\n")
cat(paste0(">>> 总细胞数: ", total_cells, "\n"))
cat(paste0(">>> T 细胞总数: ", tcell_count, "\n"))
cat(paste0(">>> T 细胞总占比: ", round(total_pct, 2), "%\n"))
cat("--------------------------------------------------\n")

# 2. 计算并绘制：基于 Tumor Stage 分组的 T 细胞占比
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
  labs(title = "T Cell Proportion across Tumor Stages", x = "Tumor Stage", y = "T Cell Percentage (%)") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
ggsave("p1_T_cell_Proportion_by_Stage.png", p1, width = 8, height = 6, dpi = 300)

# 3. 计算并绘制：样本水平 T 细胞占比与 T 细胞 Score 的相关性
sample_stats_score <- meta_data %>%
  filter(!is.na(!!sym(stage_col))) %>%
  group_by(!!sym(sample_col), !!sym(stage_col)) %>%
  summarise(
    Tcell_Pct = (sum(cluster1 == tcell_label_name, na.rm = TRUE) / n()) * 100,
    Avg_Score = mean(!!sym(score_col), na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(!is.na(Avg_Score))

p2 <- ggplot(sample_stats_score, aes(x = Avg_Score, y = Tcell_Pct)) +
  geom_point(aes(color = !!sym(stage_col)), size = 4) +
  geom_smooth(method = "lm", color = "black", fill = "lightgray") +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") + 
  theme_bw() +
  labs(title = paste0("Correlation: T Cell % vs. T-Cell ", score_col),
       x = "Average T-Cell Score (per Sample)", y = "T Cell Percentage per Sample (%)")

print(p2)
ggsave("p2_T_cell_Score_vs_Percentage.png", p2, width = 8, height = 6, dpi = 300)

# 4. 【核心交接数据】生成每个样本的 T 细胞占比表，供脚本二的癌细胞分析使用
sample_stats <- meta_data %>%
  filter(!is.na(!!sym(stage_col))) %>%
  group_by(!!sym(sample_col), !!sym(stage_col)) %>%
  summarise(
    Total_Cells = n(),                                     
    Tcell_Count = sum(cluster1 == tcell_label_name, na.rm = TRUE), 
    Tcell_Pct = (Tcell_Count / Total_Cells) * 100,         
    .groups = 'drop'
  )
write.csv(sample_stats, "Sample_Level_Tcell_Proportion.csv", row.names = FALSE)
cat(">>> T细胞分析完成！交接数据已保存至: Sample_Level_Tcell_Proportion.csv\n")