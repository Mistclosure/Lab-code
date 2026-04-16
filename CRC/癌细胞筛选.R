setwd("/mnt/disk1/qiuzerui/downloads/CRC/GSE178318")
library(data.table)
library(Seurat)
library(stringr)
library(magrittr)
library(harmony) 
library(celldex)
library(qs) 

# 1. 读取数据并创建对象
data <- Read10X("GSE178318sc") 
scRNA = CreateSeuratObject(data, min.cells = 3, project = "GSE178318", min.features = 300)

# 【核心修改】：从 Barcode 中提取病人名 (COL07 等) 并赋值给 orig.ident
# 根据你提供的 "AAACCTGAGAAACCTA_COL07_CRC"，我们取下划线分隔的第 2 部分
scRNA$orig.ident <- str_split(colnames(scRNA), "_", simplify = TRUE)[, 2]

# 计算线粒体比例
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")

### SCT 优化
scRNA <- SCTransform(scRNA, 
                     method = "glmGamPoi", 
                     vst.flavor = "v2", 
                     vars.to.regress = "percent.mt", 
                     verbose = FALSE)

### PCA
scRNA <- RunPCA(scRNA, npcs=50, verbose=FALSE)

### Harmony
# 此时 orig.ident 已经包含 COL07, COLXX 等多个水平，不再报错
scRNA <- RunHarmony(scRNA, group.by.vars="orig.ident", 
                    max.iter.harmony = 20, assay.use = "SCT")

pc.num=1:30
scRNA <- RunTSNE(scRNA, reduction="harmony", dims=pc.num) %>%
  RunUMAP(reduction="harmony", dims=pc.num)

# FindNeighbors
scRNA <- FindNeighbors(scRNA, reduction="harmony", dims = pc.num)
scRNA <- FindClusters(scRNA, resolution = 1.5)

# 读取 CSV 文件并补充注释
cluster <- read.csv("crc10x_full_c295v4_submit_cluster.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
scRNA@meta.data$cluster <- cluster$clMidwayPr[match(colnames(scRNA), cluster$sampleID)]

# 使用 qsave 保存
qsave(scRNA, "scRNA.qs")
scRNA=qread('scRNA.qs')
# 1. 加载所需的包
library(copykat)
library(dplyr)

# 2. 提取全部细胞矩阵
scRNA <- JoinLayers(scRNA, assay = "RNA")
raw_counts <- GetAssayData(scRNA, assay = "RNA", layer = "counts")

# ==================== 按 orig.ident (现在是病人名) 分批运行 CopyKAT ====================

sample_ids <- unique(scRNA$orig.ident)
pred_list <- list()

for (i in seq_along(sample_ids)) {
  current_id <- sample_ids[i]
  print(paste0("========== 正在运行第 ", i, " 组 (共 ", length(sample_ids), " 组): ", current_id, " =========="))
  
  cell_idx <- which(scRNA$orig.ident == current_id)
  sub_counts <- raw_counts[, cell_idx, drop = FALSE]
  
  sub_copykat_res <- copykat(
    rawmat = as.matrix(sub_counts), 
    id.type = "S",            
    ngene.chr = 5,            
    win.size = 25,            
    KS.cut = 0.1,             
    sam.name = paste0("copykat_", current_id), 
    distance = "euclidean", 
    n.cores = 16               
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
pred=qread("copykat_merged_pred.qs")
# 4. 提取预测结果并筛选恶性细胞
malig_barcodes <- pred$cell.names[pred$copykat.pred == "aneuploid"]

# 通过 match 直接获取对应的病人名
malig_ids <- scRNA$orig.ident[match(malig_barcodes, colnames(scRNA))]

# --- 任务 A: 统计每个 ID 的恶性细胞数量 ---
malig_counts <- as.data.frame(table(malig_ids))
colnames(malig_counts) <- c("id", "number")

cli_full <- read.csv("GSE178318_Cli.csv", header = TRUE)
meat <- merge(cli_full[,"Patient.ID",drop=FALSE], malig_counts, 
              by.x = "Patient.ID", by.y = "id", all.x = TRUE)

meat$number[is.na(meat$number)] <- 0

write.table(meat, file = "Paint_Malignant cells.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- 任务 B: 提取全部恶性细胞列表 ---
malig_df <- data.frame(Barcode = malig_barcodes, row.names = malig_barcodes)

write.table(malig_df, file = "Malignant cells.txt", 
            sep = "\t", quote = FALSE, col.names = TRUE)

# ==================== 针对恶性细胞子集的重跑流程 ====================

library(limma)
Malignant = read.table("Malignant cells.txt", header=T, sep="\t", check.names=F, row.names=1)
pbmc1 = scRNA[, rownames(Malignant)] 

pbmc1 <- JoinLayers(pbmc1,assay= 'RNA')
pbmc1[["percent.mt"]] <- PercentageFeatureSet(pbmc1, pattern = "^MT-")
pbmc1 <- SCTransform(pbmc1, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)

DefaultAssay(pbmc1) <- "RNA"
pbmc1 <- NormalizeData(pbmc1)
pbmc1 <- ScaleData(pbmc1)

DefaultAssay(pbmc1) <- "SCT" 
pbmc1 <- RunPCA(pbmc1, npcs = 50, verbose = FALSE)

# 这里同样会基于修正后的 orig.ident 进行 Harmony
pbmc1 <- RunHarmony(pbmc1, group.by.vars = "orig.ident", 
                    assay.use = "SCT", max.iter.harmony = 20)

pc.num = 1:20 
pbmc1 <- RunTSNE(pbmc1, reduction = "harmony", dims = pc.num) %>%
  RunUMAP(reduction = "harmony", dims = pc.num)

pbmc1 <- FindNeighbors(pbmc1, reduction = "harmony", dims = pc.num)
pbmc1 <- FindClusters(pbmc1, resolution = 1)

p1 <- DimPlot(pbmc1, reduction = "umap", label = TRUE) + NoLegend()
p2 <- DimPlot(pbmc1, reduction = "umap", group.by = "orig.ident")

qsave(pbmc1, "Malignant.qs")