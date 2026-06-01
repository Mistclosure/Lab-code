# coldmouse_introns 安装指南

## 环境要求

- R >= 4.3
- conda 环境: cmintro

## 创建 conda 环境

```bash
conda create -n cmintro -c conda-forge -c bioconda r-base=4.3 r-seurat r-tidyverse r-qs r-patchwork r-hdf5r r-reticulate
conda activate cmintro
```

## R 包安装

### CRAN/conda 包

```r
install.packages(c(
  "Seurat", "tidyverse", "Matrix", "qs", "patchwork",
  "ggplot2", "ggrepel", "harmony", "leidenbase", "msigdbr",
  "yaml", "hdf5r", "reticulate"
))
```

### Bioconductor 包

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment", "scDblFinder", "SingleR", "celldex",
  "clusterProfiler", "org.Mm.eg.db", "enrichplot"
))
```

### monocle2 安装

monocle2 需要特殊安装，不能通过标准 BiocManager 安装：

```r
# 安装依赖
install.packages(c("Hmisc", "igraph", "DDRTree", "densityClust"))

# 从 GitHub 安装 monocle2
devtools::install_github("cole-trapnell-lab/monocle-release")
```

## 验证安装

运行 `env/install_packages.R` 检查所有包是否安装成功。

## 包策略

- CRAN/conda: Seurat、tidyverse、Matrix、qs、patchwork、ggplot2、ggrepel、harmony、leidenbase、msigdbr、yaml、hdf5r
- Bioconductor: SingleCellExperiment、scDblFinder、SingleR、celldex、clusterProfiler、org.Mm.eg.db、enrichplot、monocle
- monocle 指 monocle2 包，不是 monocle3
