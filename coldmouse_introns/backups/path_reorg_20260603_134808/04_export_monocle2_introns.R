# ==============================================================================
# 04_export_monocle2_introns.R
# coldmouse_introns：将 Seurat 对象转换为 Monocle 2 CellDataSet
# ==============================================================================
# 输入：coldmouse_introns_sc_by_tissue_no_harmony_louvain.qs
# 输出：Monocle 2 cds 对象
# 参考：docs/original_reference/Pseudotime Trajectory Analysis_*.r
# 注意：monocle2 包不是 monocle3，转换方式不同
# ==============================================================================

# 加载配置和工具
PROJECT_HOME <- Sys.getenv("COLDMOUSE_INTRONS_HOME", unset = "/home/zerui/code/coldmouse_introns")
source(file.path(PROJECT_HOME, "R", "utils_paths.R"))
source(file.path(PROJECT_HOME, "R", "utils_io.R"))
source(file.path(PROJECT_HOME, "R", "utils_monocle2.R"))

config <- load_config()
setup_performance(config)

# 加载 R 包
library(Seurat)
library(qs)
library(monocle)
library(ggplot2)
library(dplyr)

# 设置输出目录
objects_dir <- get_results_dir(config, "objects")
files_dir <- get_results_dir(config, "files")
plots_dir <- get_results_dir(config, "plots")

# 聚类参数
cluster_method <- config$clustering$method

# 读取按组织对象（位于实验目录根层）
sc_by_tissue_file <- file.path(config$paths$output_dir,
                               paste0("coldmouse_introns_sc_by_tissue_no_harmony_", cluster_method, ".qs"))
if (!file.exists(sc_by_tissue_file)) {
  stop("未找到按组织对象文件：", sc_by_tissue_file)
}
sc_by_tissue <- load_seurat(sc_by_tissue_file)

# ------------------------------------------------------------------------------
# 1. 转换每个组织为 Monocle 2 对象
# ------------------------------------------------------------------------------
print("步骤1/2：将 Seurat 对象转换为 Monocle 2 CellDataSet...")

cds_list <- list()
for (tissue_name in names(sc_by_tissue)) {
  print(paste(">>> 转换组织:", tissue_name))
  obj <- sc_by_tissue[[tissue_name]]

  tryCatch({
    cds <- seurat_to_monocle2(obj)
    cds_list[[tissue_name]] <- cds
    print(paste("   完成:", tissue_name, "| CellDataSet 细胞数:", ncol(cds)))
  }, error = function(e) {
    message(paste("   转换组织", tissue_name, "时出错:", e$message))
  })
}

if (length(cds_list) == 0) {
  stop("没有成功转换的组织对象。")
}

# ------------------------------------------------------------------------------
# 2. 保存 Monocle 2 对象
# ------------------------------------------------------------------------------
print("步骤2/2：保存 Monocle 2 CellDataSet...")

for (tissue_name in names(cds_list)) {
  output_file <- file.path(objects_dir,
                           paste0("coldmouse_introns_monocle2_", cluster_method, "_", tissue_name, ".qs"))
  save_seurat(cds_list[[tissue_name]], output_file)
}

# 保存合并的列表
save_seurat(cds_list,
            file.path(objects_dir,
                      paste0("coldmouse_introns_monocle2_", cluster_method, "_all_tissues.qs")))

print("步骤2/2：脚本结束。Monocle 2 CellDataSet 已保存到 results/objects/。")
