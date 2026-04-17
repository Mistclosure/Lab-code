# ==============================================================================
# _targets.R
# targets pipeline 主配置文件
# ==============================================================================
library(targets)
library(tarchetypes)

# 设定项目根目录
setwd('/mnt/disk1/qiuzerui/expriments/coldmouse')

# 设置 targets 选项和环境变量
tar_option_set(
  packages = c(
    "Seurat", "tidyverse", "Matrix", "qs", "patchwork",
    "leidenbase", "scDblFinder", "SingleCellExperiment",
    "harmony", "SingleR", "celldex", "clusterProfiler",
    "org.Mm.eg.db", "enrichplot", "reticulate"
  ),
  format = "qs" # 使用 qs 处理所有内存中流转的单细胞对象，速度极快
)

# 引入我们在 R/functions.R 中定义的函数
source("R/functions.R")

# 定义数据流和依赖管线
list(
  # 追踪原始数据文件夹状态 (若文件夹内文件变化将触发下游重跑)
  tar_target(
    raw_data_dir,
    "/mnt/disk1/qiuzerui/expriments/coldmouse/sccoldmouse",
    format = "file"
  ),
  
  # Script 1: 数据加载与质控
  tar_target(
    sc_combined_qc,
    run_data_loading_and_qc(raw_data_dir)
  ),
  
  # Script 2: 聚类与全局注释
  tar_target(
    sc_by_tissue,
    run_clustering_and_global_anno(sc_combined_qc)
  ),
  
  # Script 3: 髓系细分
  tar_target(
    sc_myeloid_annotated,
    run_myeloid_clustering_and_anno(sc_by_tissue)
  ),
  
  # Script 4: Fkbp5 靶向绘图
  tar_target(
    fkbp5_targeted_res,
    run_fkbp5_targeted_analysis(sc_by_tissue)
  ),
  
  # Script 5: PBMC Cold DEA 富集
  tar_target(
    pbmc_dea_res,
    run_pbmc_dea_enrichment(sc_by_tissue)
  )
)