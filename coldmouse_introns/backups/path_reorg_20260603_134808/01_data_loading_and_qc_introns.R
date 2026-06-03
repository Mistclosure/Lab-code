# ==============================================================================
# 01_data_loading_and_qc_introns.R
# coldmouse_introns：数据读取、质量控制、核糖体基因移除和双细胞过滤
# ==============================================================================
# 输入：/mnt/disk1/qiuzerui/expriments/coldmouse_introns/rawdata/*.loom
# 输出：coldmouse_introns_sc_combined_qc_Cleaned.qs、QC summary csv、QC 图片
# 参考：docs/original_reference/data_loading_and_qc.r
# 主要修改：
#   - 读取 introns loom 文件而非 10X 矩阵
#   - 路径和参数从 config.yaml 获取
#   - QC 阈值可配置
# ==============================================================================

# 加载配置和工具
PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))
source(file.path(PROJECT_HOME, "R", "utils_qc.R"))

config <- load_config()
setup_performance(config)

# 加载 R 包
library(Seurat)
library(tidyverse)
library(Matrix)
library(qs)
library(patchwork)
library(hdf5r)
library(scDblFinder)
library(SingleCellExperiment)

# 设置输出目录
output_dir <- get_output_dir(config)
objects_dir <- get_results_dir(config, "objects")
files_dir <- get_results_dir(config, "files")
plots_dir <- get_results_dir(config, "plots")

# 读取样本映射
sample_map <- read_sample_map(config$sample_map)
sample_ids <- sample_map$sample_id

# QC 参数
qc_config <- config$qc
min_features <- qc_config$min_features
max_features <- qc_config$max_features
max_count <- qc_config$max_count
percent_mt <- qc_config$percent_mt_threshold

# 核糖体基因处理配置
ribo_config <- config$ribo
remove_ribo <- ribo_config$remove_ribo_genes
ribo_pattern <- ribo_config$ribo_pattern

# ------------------------------------------------------------------------------
# 1. 读取 introns loom 文件并创建 Seurat 对象
# ------------------------------------------------------------------------------
print("步骤1/4：读取 introns loom 矩阵并创建 Seurat 对象...")

sc_list <- list()
data_dir <- config$paths$rawdata_dir

for (i in 1:nrow(sample_map)) {
  sample <- sample_map$sample_id[i]
  loom_file <- file.path(data_dir, sample_map$loom_file[i])

  tryCatch({
    counts <- read_introns_loom(loom_file, sample = sample)
    sc_obj <- CreateSeuratObject(counts = counts, project = sample,
                                 min.cells = qc_config$min_cells,
                                 min.features = min_features)

    if (ncol(sc_obj) == 0) next

    # 添加 metadata
    sc_obj$Tissue <- sample_map$tissue[i]
    sc_obj$Group <- sample_map$group[i]
    sc_obj$Condition <- sample_map$condition[i]
    sc_obj <- calculate_qc_metrics(sc_obj, ribo_pattern = ribo_pattern)

    sc_list[[sample]] <- sc_obj
    print(paste("   样本读取完成:", sample, "| 初始细胞数:", ncol(sc_obj)))

    rm(counts, sc_obj)
    gc(verbose = FALSE)
  }, error = function(e) {
    message(paste("   处理样本时出错", sample, ":", e$message))
  })
}

if (length(sc_list) == 0) {
  stop("没有样本成功读取，请检查 loom 输入路径。")
}

sc_combined <- merge(sc_list[[1]], y = sc_list[2:length(sc_list)], add.cell.ids = names(sc_list))
sc_combined <- JoinLayers(sc_combined)

# ------------------------------------------------------------------------------
# 2. 综合质控：数值过滤 + 核糖体基因过滤
# ------------------------------------------------------------------------------
print("步骤2/4：执行 QC 过滤和核糖体基因移除...")

# 保存 QC 分位数表
qc_quantiles <- calculate_qc_quantiles(sc_combined)
write.csv(qc_quantiles, file.path(files_dir, "coldmouse_introns_qc_quantiles_before_filter.csv"),
          row.names = FALSE)

# 应用 QC 过滤
sc_combined <- apply_qc_filter(sc_combined,
                               min_features = min_features,
                               max_features = max_features,
                               max_count = max_count,
                               percent_mt = percent_mt)

print(paste("   QC 过滤后细胞数:", ncol(sc_combined)))

# 核糖体基因处理
if (remove_ribo) {
  sc_combined[["percent.rb"]] <- PercentageFeatureSet(sc_combined, pattern = ribo_pattern)
  sc_combined <- remove_ribo_genes(sc_combined, ribo_pattern = ribo_pattern)
  print(paste("   移除核糖体基因后基因数:", nrow(sc_combined)))
}

# 保存 QC 分位数表 (过滤后)
qc_quantiles_after <- calculate_qc_quantiles(sc_combined)
write.csv(qc_quantiles_after, file.path(files_dir, "coldmouse_introns_qc_quantiles_after_filter.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. 按组织拆分并分别运行 scDblFinder
# ------------------------------------------------------------------------------
print("步骤3/4：按组织拆分并运行 scDblFinder...")

sc_by_tissue <- SplitObject(sc_combined, split.by = "Tissue")

sc_by_tissue_clean <- lapply(names(sc_by_tissue), function(x) {
  run_scDblFinder(sc_by_tissue[[x]],
                  resolution = config$doublet$scDblFinder_resolution,
                  pca_dims = config$doublet$scDblFinder_pca_dims)
})
names(sc_by_tissue_clean) <- names(sc_by_tissue)
sc_by_tissue_clean <- sc_by_tissue_clean[!sapply(sc_by_tissue_clean, is.null)]

if (length(sc_by_tissue_clean) == 0) {
  stop("QC/去双胞后没有剩余组织对象。")
}

# ------------------------------------------------------------------------------
# 4. 合并、生成 QC 图片并保存
# ------------------------------------------------------------------------------
print("步骤4/4：合并清洗后的对象并保存 introns QC 输出...")

if (length(sc_by_tissue_clean) == 1) {
  sc_final <- sc_by_tissue_clean[[1]]
} else {
  sc_final <- merge(sc_by_tissue_clean[[1]], y = sc_by_tissue_clean[2:length(sc_by_tissue_clean)])
}
sc_final <- JoinLayers(sc_final)

# 生成 QC 图片
generate_qc_plots(sc_final, output_dir = plots_dir, prefix = "coldmouse_introns_sc_combined_qc")

# 保存 QC summary
qc_summary <- data.frame(
  sample = names(sc_list),
  initial_cells = sapply(sc_list, ncol),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(files_dir, "coldmouse_introns_qc_summary.csv"), row.names = FALSE)

# 保存最终对象
output_file <- file.path(objects_dir, "coldmouse_introns_sc_combined_qc_Cleaned.qs")
save_seurat(sc_final, output_file)

print(paste("完成。最终高质量细胞数:", ncol(sc_final)))
print(paste("已保存:", output_file))
