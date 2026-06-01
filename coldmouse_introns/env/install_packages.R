# ==============================================================================
# install_packages.R: coldmouse_introns 项目 R 包安装脚本
# ==============================================================================
# 此脚本用于安装项目所需的所有 R 包。
# 不要在生产脚本中自动安装，所有安装操作通过此脚本完成。
# ==============================================================================

# CRAN/conda 包
cran_packages <- c(
  "Seurat",
  "tidyverse",
  "Matrix",
  "qs",
  "patchwork",
  "ggplot2",
  "ggrepel",
  "harmony",
  "leidenbase",
  "msigdbr",
  "yaml",
  "hdf5r",
  "reticulate"
)

# Bioconductor 包
bioc_packages <- c(
  "SingleCellExperiment",
  "scDblFinder",
  "SingleR",
  "celldex",
  "clusterProfiler",
  "org.Mm.eg.db",
  "enrichplot"
)

# 安装 CRAN 包
message("安装 CRAN 包...")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("安装: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    message("已安装: ", pkg)
  }
}

# 安装 Bioconductor 包
message("安装 Bioconductor 包...")
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
for (pkg in bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message("安装: ", pkg)
    BiocManager::install(pkg, ask = FALSE)
  } else {
    message("已安装: ", pkg)
  }
}

# 检查 monocle2 (需要特殊处理)
message("检查 monocle2...")
if (!require("monocle", quietly = TRUE)) {
  message("monocle2 未安装。请参考 docs/install.md 手动安装。")
} else {
  message("已安装: monocle2")
}

message("包安装检查完成。")
