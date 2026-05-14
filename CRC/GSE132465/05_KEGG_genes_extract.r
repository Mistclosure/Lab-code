# ==============================================================================
# 提取特定通路及其对应的基因列表 (Long Format)
# ==============================================================================

library(dplyr)
library(tidyr)

# 1. 设置工作目录 (请确保与你之前的脚本一致)
WORK_DIR <- '/mnt/disk1/qiuzerui/downloads/CRC/GSE132465'
FILES_DIR <- "files"
setwd(WORK_DIR)

# 2. 读取之前生成的 compareCluster 结果文件
input_file <- file.path(FILES_DIR, "Target_Modules_compareCluster_KEGG_signed.csv")

if (!file.exists(input_file)) {
  stop("找不到输入文件，请确认路径是否正确或是否已运行上一个分析脚本。")
}

kegg_data <- read.csv(input_file, stringsAsFactors = FALSE)

# 3. 定义你想要筛选的四个目标通路
target_pathways <- c(
  "Cell cycle", 
  "Rap1 signaling pathway", 
  "Proteoglycans in cancer", 
  "Gap junction"
)

# 4. 筛选数据并进行“一变多”转换
# 逻辑：筛选指定通路 -> 只保留 Description (通路名) 和 geneID (基因列表) -> 将 geneID 按 "/" 拆分
pathway_gene_map <- kegg_data %>%
  filter(Description %in% target_pathways) %>%
  select(Description, geneID) %>%
  # 将 "Gene1/Gene2/Gene3" 转换为多行，每行一个基因
  separate_rows(geneID, sep = "/") %>%
  # 去重（防止同一个基因在同一个通路中因为不同模块重复出现）
  distinct() %>%
  # 按照要求重命名列名并调整顺序
  rename(Pathway = Description, Gene = geneID) %>%
  select(Gene, Pathway)

# 5. 输出结果
output_file <- file.path(FILES_DIR, "CRC_Proliferation_Invasion_Metastasis_Genes.csv")
write.csv(pathway_gene_map, file = output_file, row.names = FALSE, quote = FALSE)

cat("\n================================================================\n")
cat("任务完成！\n")
cat("筛选通路:", paste(target_pathways, collapse = ", "), "\n")
cat("结果已保存至:", output_file, "\n")
cat("总计提取条目数:", nrow(pathway_gene_map), "\n")
cat("================================================================\n")