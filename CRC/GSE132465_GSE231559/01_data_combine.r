# ==============================================================================
# 单细胞数据合并与 Harmony 重新整合流程
# 数据集: GSE132465, GSE231559
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
# 添加 Dataset 标签
sc132465$Dataset <- "GSE132465"
sc231559$Dataset <- "GSE231559"

# 合并对象并为 cell barcode 添加前缀以防重名
merged_sc <- merge(x = sc132465, 
                   y = sc231559,  
                   project = "CRC_Merged")

# 清理内存
rm(sc132465, sc231559)
gc()

# 4. 全局标准化与降维预处理 ----------------------------------------------------
cat("正在进行标准化与 PCA 降维...\n")
DefaultAssay(merged_sc) <- "RNA"

# 对数据进行对数标准化，设置缩放因子为 1000000
merged_sc <- NormalizeData(merged_sc, normalization.method = "LogNormalize", scale.factor = 1000000)

# 寻找 2000 个高可变基因
merged_sc <- FindVariableFeatures(merged_sc, selection.method = "vst", nfeatures = 2000)

cat("正在进行 ScaleData (这可能需要一些时间)...\n")
# 【修正注释以匹配功能】：未指定 features 参数，Seurat 默认仅对上述 2000 个高可变基因进行缩放 (Scale)
merged_sc <- ScaleData(merged_sc, verbose = FALSE)

# 使用高可变基因进行主成分分析 (PCA)
merged_sc <- RunPCA(merged_sc, npcs = 30, verbose = FALSE)

# 5. 运行 Harmony 消除批次效应 -------------------------------------------------
cat("正在运行 Harmony 进行批次校正...\n")
# 【修正注释以匹配功能】：group.by.vars = "orig.ident" 表明实际上是在消除不同样本(Patient/Sample)之间的批次效应，而不仅仅是数据集间的批次效应
merged_sc <- RunHarmony(merged_sc, 
                        group.by.vars = "orig.ident", 
                        plot_convergence = TRUE)
                        # ==============================================================================


# 6. 使用校正后的数据进行 UMAP 降维与细胞聚类 ----------------------------------
cat("正在运行 UMAP 与聚类...\n")
# 基于 Harmony 降维结果计算 UMAP 和构建邻接图
merged_sc <- RunUMAP(merged_sc, reduction = "harmony", dims = 1:30)
merged_sc <- FindNeighbors(merged_sc, reduction = "harmony", dims = 1:30)
merged_sc <- FindClusters(merged_sc, resolution = 0.5) 

# 7. 合并临床信息 (安全模式) ---------------------------------------------------
cat("正在合并临床信息...\n")
# 处理特定的样本命名问题
merged_sc$orig.ident[merged_sc$orig.ident %in% c('L8T1','L8T2')] <- 'L8T'

# 读取临床信息表
cli <- read.csv('GSE132465_GSE231559_Cli.csv', header = TRUE, check.names = FALSE)

# 安全检查，防止行数因合并增加
if (any(duplicated(cli$Patient))) {
  stop("临床信息表中的 Patient 列有重复项，请先去重后再运行！")
}

# 提取 meta.data 并保留行名为单独一列
meta_df <- merged_sc@meta.data
meta_df$Cell_Barcode <- rownames(meta_df)

# 使用 left_join 按照 Patient 列合并临床信息
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

print(p1 + p2) 

# 9. 保存全局整合结果并提取癌细胞 ----------------------------------------------
cat("正在保存合并并整合后的 Seurat 对象...\n")
qsave(merged_sc, "Merged_scRNA.qs")
merged_sc = qread('Merged_scRNA.qs')

cat("正在提取恶性细胞子集...\n")
Malignant132 <- read.table("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/Malignant cells.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
Malignant231 <- read.table("/mnt/disk1/qiuzerui/downloads/CRC/GSE231559/Malignant cells.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# 根据提供的 TXT 提取恶性细胞
malig <- merged_sc[, colnames(merged_sc) %in% c(rownames(Malignant132),rownames(Malignant231))] 


# ==============================================================================
# 10. 恶性细胞子集重跑 RNA 流程 (参照全局流程参数)
# ==============================================================================
cat("正在对恶性细胞子集重新进行预处理与降维...\n")

# 对恶性细胞子集重新进行对数标准化，比例因子保持一致 (1000000)
malig <- NormalizeData(malig, normalization.method = "LogNormalize", scale.factor = 1000000)

# 重新寻找仅针对恶性细胞内部变异的 2000 个高可变基因
malig <- FindVariableFeatures(malig, selection.method = "vst", nfeatures = 2000)

# 对恶性细胞子集的高可变基因进行缩放
malig <- ScaleData(malig, verbose = FALSE)

# 针对恶性细胞子集重新进行 PCA 主成分分析
malig <- RunPCA(malig, npcs = 30, verbose = FALSE)

# 重新运行 Harmony，消除恶性细胞子集在不同样本 (orig.ident) 间的批次效应
cat("正在对恶性细胞子集运行 Harmony...\n")
malig <- RunHarmony(malig, 
                    group.by.vars = "orig.ident", 
                    plot_convergence = TRUE)

# 对恶性细胞子集重新进行 UMAP 降维与细胞分群
cat("正在对恶性细胞子集运行 UMAP 与聚类...\n")
malig <- RunUMAP(malig, reduction = "harmony", dims = 1:30)
malig <- FindNeighbors(malig, reduction = "harmony", dims = 1:30)
malig <- FindClusters(malig, resolution = 0.5)

# 保存最终的癌细胞特异性 Seurat 对象
cat("流程运行完毕，正在保存恶性细胞 Seurat 对象...\n")
qsave(malig, 'Merged_malig.qs')