# ==============================================================================
# 单细胞数据合并与 Harmony 重新整合流程
# 数据集: GSE132465, GSE231559
# 更新: Scale全部基因 & 安全合并临床信息
# ==============================================================================

# 1. 载入所需的 R 包 -----------------------------------------------------------
library(Seurat)
library(qs)
library(harmony)
library(dplyr)
library(ggplot2)

# 2. 设置工作目录并读取数据 ----------------------------------------------------
setwd('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465_GSE231559')

cat("正在读取数据...\n")
sc132465 <- qread('/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/scRNA_RNA_assay.qs')
sc231559 <- qread('/mnt/disk1/qiuzerui/downloads/CRC/GSE231559/scRNA_RNA_assay.qs')

# 3. 添加数据集标签并合并对象 --------------------------------------------------
cat("正在合并 Seurat 对象...\n")
# 为后续 Harmony 批次校正添加 Dataset 标签
sc132465$Dataset <- "GSE132465"
sc231559$Dataset <- "GSE231559"

# 合并对象并为 cell barcode 添加前缀以防重名
merged_sc <- merge(x = sc132465, 
                   y = sc231559, 
                   add.cell.ids = c("GSE132465", "GSE231559"), 
                   project = "CRC_Merged")

# 清理内存
rm(sc132465, sc231559)
gc()

# 4. 重新进行全局标准化与降维预处理 --------------------------------------------
cat("正在进行标准化与 PCA 降维...\n")
DefaultAssay(merged_sc) <- "RNA"

merged_sc <- NormalizeData(merged_sc, normalization.method = "LogNormalize", scale.factor = 1000000)
merged_sc <- FindVariableFeatures(merged_sc, selection.method = "vst", nfeatures = 2000)

cat("正在对所有基因进行 ScaleData (这可能需要一些时间)...\n")
# 【修改点 1】：ScaleData 改为全部基因，为后续 CRDscore 打分做准备
merged_sc <- ScaleData(merged_sc, verbose = FALSE)

merged_sc <- RunPCA(merged_sc, npcs = 30, verbose = FALSE)

# 5. 运行 Harmony 消除数据集间的批次效应 ---------------------------------------
cat("正在运行 Harmony 进行批次校正...\n")
merged_sc <- RunHarmony(merged_sc, 
                        group.by.vars = "orig.ident", 
                        plot_convergence = TRUE)

# 6. 使用校正后的数据进行 UMAP 降维与细胞聚类 ----------------------------------
cat("正在运行 UMAP 与聚类...\n")
merged_sc <- RunUMAP(merged_sc, reduction = "harmony", dims = 1:30)
merged_sc <- FindNeighbors(merged_sc, reduction = "harmony", dims = 1:30)
merged_sc <- FindClusters(merged_sc, resolution = 0.5) 

# 7. 合并临床信息 (安全模式) ---------------------------------------------------
cat("正在合并临床信息...\n")
# 处理特定的样本命名问题
merged_sc$orig.ident[merged_sc$orig.ident %in% c('L8T1','L8T2')] <- 'L8T'

# 读取临床信息表
cli <- read.csv('GSE132465_GSE231559_Cli.csv', header = TRUE, check.names = FALSE)

# 【修改点 2】：安全合并临床信息，防止打乱细胞顺序
if (any(duplicated(cli$Patient))) {
  stop("临床信息表中的 Patient 列有重复项，请先去重后再运行！")
}

# 提取 meta.data 并保留行名为单独一列
meta_df <- merged_sc@meta.data
meta_df$Cell_Barcode <- rownames(meta_df)

# 使用 left_join 合并临床信息
meta_merged <- left_join(meta_df, cli, by = c("orig.ident" = "Patient"))

# 严格恢复原来的行名和细胞顺序
rownames(meta_merged) <- meta_merged$Cell_Barcode
meta_merged <- meta_merged[rownames(merged_sc@meta.data), ]

# 去除临时的 Barcode 列并将数据覆盖回 Seurat 对象
meta_merged$Cell_Barcode <- NULL
merged_sc@meta.data <- meta_merged

# 8. 结果可视化与检查 ----------------------------------------------------------
cat("正在生成 UMAP 图像...\n")
p1 <- DimPlot(merged_sc, reduction = "umap", group.by = "Dataset") + ggtitle("Batch Effect check (By Dataset)")
p2 <- DimPlot(merged_sc, reduction = "umap", group.by = "seurat_clusters", label = TRUE) + ggtitle("By Clusters")

# 显示图片 (如果在非图形界面运行，请注释掉 print，改用 ggsave)
print(p1 + p2) 
# ggsave("UMAP_Harmony_Check.pdf", plot = p1 + p2, width = 12, height = 5)

# 9. 保存最终结果 --------------------------------------------------------------
cat("正在保存合并并整合后的 Seurat 对象...\n")
qsave(merged_sc, "Merged_scRNA.qs")

cat("流程运行完毕！\n")