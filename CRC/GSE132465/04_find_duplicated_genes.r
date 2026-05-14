# ==========================================
# 1. 读取两个 CSV 文件
# ==========================================
# 请将引号内的文件名替换为你实际的文件路径（例如 "C:/data/geneset1.csv"）
data1 <- read.csv("/mnt/disk1/qiuzerui/downloads/CRC/signature/ciliopathy_genes.csv", stringsAsFactors = FALSE)
data2 <- read.csv("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/target_gene_df.csv", stringsAsFactors = FALSE)
#data2 <- read.table("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/files/Target_Modules_Signed_单独列表/Module_blue_signed_genes.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# ==========================================
# 2. 提取基因名称列
# ==========================================
# 假设你的基因名称所在的列名叫 "GeneSymbol"，请根据你的 CSV 实际列名进行修改。
# 如果你的基因名在第一列，也可以直接用 genes1 <- data1[, 1]
genes1 <- data1[, 1]
genes2 <- data2[, 1]

# 清理可能存在的空白字符或去重（良好习惯）
genes1 <- unique(trimws(genes1))
genes2 <- unique(trimws(genes2))

# ==========================================
# 3. 获取重合基因（取交集）
# ==========================================
overlapping_genes <- intersect(genes1, genes2)

# 在控制台打印重合基因的数量
cat("第一个基因集数量:", length(genes1), "\n")
cat("第二个基因集数量:", length(genes2), "\n")
cat("重合的基因数量:", length(overlapping_genes), "\n")

# ==========================================
# 4. 将重合基因保存为新的 CSV 文件
# ==========================================
# 将结果转换为数据框格式以便导出
result_df <- data.frame(Overlapping_Genes = overlapping_genes)

# 导出为 CSV（row.names = FALSE 表示不保存行号）
write.csv(result_df, "overlapping_genes_result.csv", row.names = FALSE)
cat("重合基因已成功保存至当前工作目录下的 'overlapping_genes_result.csv'\n")

# ==========================================
# 5. 可选进阶：绘制韦恩图 (Venn Diagram) 直观展示
# ==========================================
# 如果你还没有安装这两个包，请取消下面这行代码的注释来安装：
# install.packages(c("ggVennDiagram", "ggplot2"))

library(ggVennDiagram)
library(ggplot2)

# 将两个基因集放入一个列表中
gene_list <- list(
  Set1 = genes1,  # 可以将 "Set1" 改为你的数据集名称，如 "Module_Brown"
  Set2 = genes2
)

# 绘图
p <- ggVennDiagram(gene_list) + 
  scale_fill_gradient(low = "white", high = "red") +  # 设置颜色渐变
  labs(title = "Overlap of Two Gene Sets") +          # 设置标题
  theme(plot.title = element_text(hjust = 0.5))       # 标题居中

# 显示图形
print(p)
data3 = data2[data2$Gene %in% result_df[,1],]
