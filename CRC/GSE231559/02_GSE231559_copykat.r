# ==============================================================================
# 脚本2：标准 RNA 流程降维聚类 + CopyKAT 恶性细胞鉴定 + 恶性细胞亚群重分析
# ==============================================================================

setwd("/mnt/disk1/qiuzerui/downloads/CRC/GSE231559")
library(data.table)
library(Seurat)
library(stringr)
library(magrittr)
library(harmony) 
library(celldex)
library(qs) 
library(copykat)
library(dplyr)
library(limma)

# ------------------------------------------------------------------------------
# 1. 读取脚本 1 预处理好的数据
# ------------------------------------------------------------------------------
print("🚀 正在读取脚本 1 生成的质控后数据...")
# 读取上一步保存的 Cleaned 文件
scRNA <- qread("sc_combined_qc_Cleaned.qs")

# 优化：Seurat V5 merge 后必须 JoinLayers 才能进行后续标准归一化
scRNA <- subset(scRNA, subset = Type == "Tumor")
scRNA <- JoinLayers(scRNA)
# 计算线粒体比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

# ------------------------------------------------------------------------------
# 2. 标准 RNA 分析流程 (替换原有的 SCT)
# ------------------------------------------------------------------------------
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
# 优化：在 ScaleData 中回归线粒体基因，保留原 SCT 逻辑
scRNA <- ScaleData(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)

### PCA
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

### Harmony (适配 RNA assay)
scRNA <- RunHarmony(scRNA, group.by.vars = "orig.ident", 
                    max.iter.harmony = 20, assay.use = "RNA")

pc.num = 1:30
scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = pc.num) %>%
  RunUMAP(reduction = "harmony", dims = pc.num)

# FindNeighbors
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 1.5)

# 使用 qsave 保存整体对象
qsave(scRNA, "scRNA_RNA_assay.qs")

# ------------------------------------------------------------------------------
# 3. 按 orig.ident 分批运行 CopyKAT
# ------------------------------------------------------------------------------
# 提取全部细胞的原始矩阵用于 CopyKAT
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

sample_ids <- unique(scRNA$orig.ident)
pred_list <- list()

for (i in seq_along(sample_ids)) {
  current_id <- sample_ids[i]
  print(paste0("========== 正在运行第 ", i, " 组 (共 ", length(sample_ids), " 组): ", current_id, " =========="))
  
  cell_idx <- which(scRNA$orig.ident == current_id)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]
  
  # 优化：跳过细胞数过少的样本，防止报错中断
  if(ncol(sub_counts) < 50) {
    print(paste0("警告：样本 ", current_id, " 细胞数过少 (", ncol(sub_counts), ")，跳过 CopyKAT。"))
    next
  }
  
  sub_copykat_res <- copykat(
    rawmat = as.matrix(sub_counts), 
    id.type = "S",            
    ngene.chr = 5,            
    win.size = 25,            
    KS.cut = 0.1,             
    sam.name = paste0("copykat_", current_id), 
    distance = "euclidean", 
    n.cores = 24               
  )
  
  if (!is.null(sub_copykat_res) && "prediction" %in% names(sub_copykat_res)) {
    pred_list[[current_id]] <- as.data.frame(sub_copykat_res$prediction)
  } else {
    print(paste0("警告：样本 ", current_id, " 未能返回有效预测结果，已跳过。"))
  }
  
  rm(sub_counts, sub_copykat_res)
  gc()
}

print("========== 所有分组运行完毕，正在合并预测结果 ==========")

pred <- do.call(rbind, pred_list)
qsave(pred, "copykat_merged_pred.qs")

# ------------------------------------------------------------------------------
# 4. 提取预测结果并匹配临床信息
# ------------------------------------------------------------------------------
# 筛选恶性细胞 (aneuploid)
malig_barcodes <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 通过 match 直接获取对应的样本名
malig_ids <- scRNA$orig.ident[match(malig_barcodes, colnames(scRNA))]

# --- 任务 A: 统计每个 ID 的恶性细胞数量 ---
malig_counts <- as.data.frame(table(malig_ids))
colnames(malig_counts) <- c("id", "number")

cli_full <- read.csv("GSE178341_crc10x_full_c295v4_submit_metatables.csv", header = TRUE)

meat <- merge(cli_full[, c("ID", "Original_Sample_ID")], malig_counts, 
              by.x = "Original_Sample_ID", by.y = "id", all.x = TRUE)

meat$number[is.na(meat$number)] <- 0
meat <- meat[, c("ID", "Original_Sample_ID", "number")]
colnames(meat)[1] <- "id"

write.table(meat, file = "Paint_Malignant cells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 任务 B: 提取全部恶性细胞列表 ---
malig_df <- data.frame(Barcode = malig_barcodes, row.names = malig_barcodes)

write.table(malig_df, file = "Malignant cells.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# 5. 针对恶性细胞子集的重跑流程 (替换原有的 SCT)
# ------------------------------------------------------------------------------
Malignant <- read.table("Malignant cells.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
pbmc1 <- scRNA[, rownames(Malignant)] 

# 确保提取后的对象图层合并
pbmc1 <- JoinLayers(pbmc1)
DefaultAssay(pbmc1) <- "RNA"

# 标准 RNA 分析流程
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
pbmc1 <- ScaleData(pbmc1, vars.to.regress = "percent.mt", verbose = FALSE)

pbmc1 <- RunPCA(pbmc1, npcs = 50, verbose = FALSE)

# 使用 RNA 层跑 Harmony
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "orig.ident", 
                    assay.use = "RNA", max.iter.harmony = 20)

pc.num = 1:20 
pbmc1 <- RunTSNE(pbmc1, reduction = "harmony", dims = pc.num) %>%
  RunUMAP(reduction = "harmony", dims = pc.num)

pbmc1 <- FindNeighbors(pbmc1, reduction = "harmony", dims = pc.num)
pbmc1 <- FindClusters(pbmc1, resolution = 1)

p1 <- DimPlot(pbmc1, reduction = "umap", label = TRUE) + NoLegend()
p2 <- DimPlot(pbmc1, reduction = "umap", group.by = "orig.ident")

qsave(pbmc1, "Malignant_RNA_assay.qs")