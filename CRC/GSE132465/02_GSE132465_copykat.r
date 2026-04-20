# ==============================================================================
# 脚本2：GSE132465 标准 RNA 流程降维聚类 + CopyKAT 恶性细胞鉴定 + 亚群重分析
# ==============================================================================

setwd("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465") # 修改为新路径
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
# 1. 读取脚本 1 预处理好的纯净数据
# ------------------------------------------------------------------------------
print("🚀 正在读取脚本 1 生成的质控后数据...")
scRNA <- qread("sc_combined_qc_Cleaned.qs")

# 【适配修改】GSE132465 区分肿瘤/正常的列名为 'Class'，提取 Tumor 细胞
if ("Class" %in% colnames(scRNA@meta.data)) {
  scRNA <- subset(scRNA, subset = Class == "Tumor")
} else {
  warning("⚠️ 未找到 'Class' 列，请确认提取 Tumor 的 metadata 字段。")
}

# 必须 JoinLayers 才能进行后续标准归一化
scRNA <- JoinLayers(scRNA)

# 计算线粒体比例
if (any(grepl("^MT-", rownames(scRNA)))) {
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
} else {
  scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")
}

# ------------------------------------------------------------------------------
# 2. 标准 RNA 分析流程 
# ------------------------------------------------------------------------------
print("🚀 正在运行标准降维与 Harmony 批次效应去除...")
DefaultAssay(scRNA) <- "RNA"
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)
scRNA <- ScaleData(scRNA, vars.to.regress = "percent.mt", verbose = FALSE)

### PCA
scRNA <- RunPCA(scRNA, npcs = 50, verbose = FALSE)

### Harmony (【适配修改】按 Patient 进行批次整合，而非 orig.ident)
scRNA <- RunHarmony(scRNA, group.by.vars = "Patient", 
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
# 3. 按 Patient 分批运行 CopyKAT
# ------------------------------------------------------------------------------
# 提取全部细胞的原始矩阵用于 CopyKAT
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

# 【适配修改】按 Patient 拆分运行 CopyKAT
patient_ids <- unique(scRNA$Patient)
pred_list <- list()

for (i in seq_along(patient_ids)) {
  current_id <- patient_ids[i]
  print(paste0("========== 正在运行 CopyKAT 第 ", i, " 组 (共 ", length(patient_ids), " 患者): ", current_id, " =========="))
  
  cell_idx <- which(scRNA$Patient == current_id)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]
  
  if(ncol(sub_counts) < 50) {
    print(paste0("   ⚠️ 警告：患者 ", current_id, " 细胞数过少 (", ncol(sub_counts), ")，跳过 CopyKAT。"))
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
    print(paste0("   ⚠️ 警告：患者 ", current_id, " 未能返回有效预测结果，已跳过。"))
  }
  
  # 【优化】极其重要的内存释放
  rm(sub_counts, sub_copykat_res)
  gc()
}

print("========== 所有患者 CopyKAT 运行完毕，正在合并预测结果 ==========")

pred <- do.call(rbind, pred_list)
qsave(pred, "copykat_merged_pred.qs")
pred = qread('copykat_merged_pred.qs')

# ------------------------------------------------------------------------------
# 4. 提取预测结果并匹配临床信息
# ------------------------------------------------------------------------------
# 筛选恶性细胞 (aneuploid)
malig_barcodes <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 【适配修改】虽然是按 Patient 跑的，但统计数量时，我们需要对应到 orig.ident 
# 因为你后续的临床表合并逻辑是匹配 orig.ident (Tumor ID)
malig_ids <- scRNA$orig.ident[match(malig_barcodes, colnames(scRNA))]

# --- 任务 A: 统计每个 orig.ident 的恶性细胞数量 ---
malig_counts <- as.data.frame(table(malig_ids))
colnames(malig_counts) <- c("id", "number")

# 【适配修改】读取 GSE132465 的临床信息
cli_full <- read.csv("GSE132465_Cli.csv", header = TRUE, check.names = FALSE)

# 【适配修改】根据你的历史代码，临床信息的匹配列为 'Tumor'
meat <- merge(cli_full, malig_counts, 
              by.x = "Tumor", by.y = "id", all.x = TRUE)

meat$number[is.na(meat$number)] <- 0
# 仅保留你需要的标识列和数量，重命名回 id (如果你的下游代码需要)
meat <- meat[, c("Tumor", "number")] 
colnames(meat)[1] <- "id"

write.table(meat, file = "Paint_Malignant cells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 任务 B: 提取全部恶性细胞列表 ---
malig_df <- data.frame(Barcode = malig_barcodes, row.names = malig_barcodes)

write.table(malig_df, file = "Malignant cells.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# 5. 针对恶性细胞子集的重跑流程 
# ------------------------------------------------------------------------------
print("🚀 正在提取恶性细胞子集并重新聚类降维...")
Malignant <- read.table("Malignant cells.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
pbmc1 <- scRNA[, rownames(Malignant)] 

pbmc1 <- JoinLayers(pbmc1)
DefaultAssay(pbmc1) <- "RNA"

# 标准 RNA 分析流程
pbmc1 <- NormalizeData(pbmc1, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc1 <- FindVariableFeatures(pbmc1, selection.method = "vst", nfeatures = 2000)
pbmc1 <- ScaleData(pbmc1, vars.to.regress = "percent.mt", verbose = FALSE)

pbmc1 <- RunPCA(pbmc1, npcs = 50, verbose = FALSE)

# 【适配修改】同样使用 Patient 作为批次效应变量
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "Patient", 
                    assay.use = "RNA", max.iter.harmony = 20)

pc.num = 1:20 
pbmc1 <- RunTSNE(pbmc1, reduction = "harmony", dims = pc.num) %>%
  RunUMAP(reduction = "harmony", dims = pc.num)

pbmc1 <- FindNeighbors(pbmc1, reduction = "harmony", dims = pc.num)
pbmc1 <- FindClusters(pbmc1, resolution = 1)

p1 <- DimPlot(pbmc1, reduction = "umap", label = TRUE) + NoLegend()
# 【适配修改】检查批次效应去除情况
p2 <- DimPlot(pbmc1, reduction = "umap", group.by = "Patient")

qsave(pbmc1, "Malignant_RNA_assay.qs")
print("🎉 恶性细胞提取与重分析完成！")