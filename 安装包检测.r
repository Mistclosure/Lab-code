# ==============================================================================
# 自动依赖管理：检测缺失包并调用 Conda 安装
# ==============================================================================

# 1. 定义需要用到的包 (格式: "R中包名" = "Conda中包名")
# 注意：Bioconductor 包通常以 bioconductor- 开头，普通包通常以 r- 开头
required_pkgs_map <- c(
  "Seurat"               = "r-seurat",
  "tidyverse"            = "r-tidyverse",
  "patchwork"            = "r-patchwork",
  "Matrix"               = "r-matrix",
  "scales"               = "r-scales",
  "openxlsx"             = "r-openxlsx",
  "ggplot2"              = "r-ggplot2",
  "pheatmap"             = "r-pheatmap",
  "scDblFinder"          = "bioconductor-scdblfinder",
  "SingleCellExperiment" = "bioconductor-singlecellexperiment",
  "BiocManager"          = "r-biocmanager"
)

# 2. 定义检测安装函数
check_and_install <- function(pkg_map) {
  # 获取当前已安装的包列表
  installed <- installed.packages()[, "Package"]
  
  # 找出缺失的包
  missing_pkgs <- names(pkg_map)[!names(pkg_map) %in% installed]
  
  if (length(missing_pkgs) > 0) {
    message("⚠️ 检测到以下包缺失，正在尝试通过 Conda 自动安装: ")
    print(missing_pkgs)
    
    # 提取对应的 Conda 包名
    conda_pkgs_to_install <- pkg_map[missing_pkgs]
    
    # 构建 Conda 安装命令 (使用 mamba 加速，如果没有 mamba 会报错，可改为 conda)
    # -y 表示自动确认，-c 指定频道
    cmd <- paste0(
      "mamba install -y -c bioconda -c conda-forge ", 
      paste(conda_pkgs_to_install, collapse = " ")
    )
    
    message("🚀 正在执行系统命令: ", cmd)
    
    # 调用系统命令执行安装
    exit_code <- system(cmd)
    
    if (exit_code == 0) {
      message("✅ 所有缺失包已通过 Conda 安装成功！")
    } else {
      stop("❌ Conda 安装失败，请查看上方报错信息，或尝试手动在终端运行该命令。")
    }
    
  } else {
    message("✅ 环境检查通过：所有必需的 R 包都已安装。")
  }
}

# 3. 执行检查
check_and_install(required_pkgs_map)

# ==============================================================================
# 下面开始正式加载包...
# ==============================================================================
library(Seurat)
library(tidyverse)
# ... (继续你的脚本)