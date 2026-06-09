# Monocle2 DDRTree 版本执行说明

## 脚本位置

DDRTree 版本脚本：

```text
/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_monocle2_ddrtree.r
```

来源脚本：

```text
/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_monocle2_umap.r
```

## 分析目标

该脚本用于 PBMC Monocyte 与 Aorta Macrophage 的合并 monocle2 拟时序分析。

与 UMAP 版本的区别：

1. UMAP 版本使用 monocle2 计算 pseudotime，但最终用 Seurat UMAP 坐标作图。
2. DDRTree 版本使用 monocle2 的 `plot_cell_trajectory()` 输出 DDRTree component 轨迹图。
3. DDRTree 版本不再运行 Harmony 和 Seurat UMAP，因为最终主图不使用 UMAP 坐标。
4. 两个版本保留相同的细胞筛选逻辑、monocle2 CellDataSet 构建逻辑、ordering genes 逻辑和 Ly6c2 root_state 逻辑。

## 输入文件

脚本工作目录固定为：

```text
/mnt/disk1/qiuzerui/expriments/coldmouse/
```

需要存在：

```text
pbmc_corrected.qs
aorta_corrected.qs
```

## 细胞筛选逻辑

PBMC：

```r
pbmc_sub <- subset(pbmc, new_clusters %in% c('3', '12'))
pbmc_sub$cell_type <- ifelse(pbmc_sub$new_clusters == '12',
                              "Monocyte", "Dendritic cells")
pbmc_sub <- subset(pbmc_sub, new_clusters == '12')
```

Aorta：

```r
aorta_target <- subset(aorta, seurat_clusters %in% c('1','6','10','14'))
aorta_target$cell_type <- aorta_target$SingleR.labels
```

合并后去掉 TN_30C：

```r
merged_obj <- subset(merged_obj, Group == 'TN_30C', invert = TRUE)
```

## monocle2 核心流程

脚本使用：

```r
newCellDataSet()
estimateSizeFactors()
estimateDispersions()
detectGenes()
setOrderingFilter()
reduceDimension(method = "DDRTree")
orderCells()
```

ordering genes 优先来自：

```r
VariableFeatures(merged_obj)
```

如果 variable features 少于 100 个，则回退到：

```r
num_cells_expressed >= 10
```

拟时序起点使用 Ly6c2 高表达细胞所在最多的 monocle2 `State`：

```r
root_state <- get_root_state_by_marker(cds, marker_gene = "Ly6c2")
cds <- orderCells(cds, root_state = root_state)
```

## 输出文件

运行后会在 `/mnt/disk1/qiuzerui/expriments/coldmouse/` 下创建或使用：

```text
files/
pictures/
```

对象输出：

```text
files/script_06_monocle2_ddrtree_standard.qs
files/script_06_monocle2_ddrtree_standard.rds
```

表格输出：

```text
files/script_06_monocle2_ddrtree_pseudotime_metadata.csv
files/script_06_monocle2_ddrtree_state_group_table.csv
files/script_06_monocle2_ddrtree_state_source_table.csv
```

图片输出：

```text
pictures/monocle2_ddrtree_source.png
pictures/monocle2_ddrtree_pseudotime.png
pictures/monocle2_ddrtree_state.png
pictures/monocle2_ddrtree_group.png
pictures/monocle2_ddrtree_cell_type.png
pictures/monocle2_ddrtree_split.png

pictures/monocle2_ddrtree_source.pdf
pictures/monocle2_ddrtree_pseudotime.pdf
pictures/monocle2_ddrtree_state.pdf
pictures/monocle2_ddrtree_group.pdf
pictures/monocle2_ddrtree_cell_type.pdf
pictures/monocle2_ddrtree_split.pdf
```

## 推荐运行环境

优先使用已有 R 环境：

```bash
conda activate cmintro
```

如果该环境不可用，再检查其他 R 环境。

## 运行前检查

```bash
df -h
free -h
nproc
which R
conda env list
```

检查输入文件：

```bash
ls -lh /mnt/disk1/qiuzerui/expriments/coldmouse/pbmc_corrected.qs
ls -lh /mnt/disk1/qiuzerui/expriments/coldmouse/aorta_corrected.qs
```

## 只做语法检查

不运行完整分析，只检查 R 语法：

```bash
/home/zerui/miniconda3/envs/cmintro/bin/Rscript -e 'parse("/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_monocle2_ddrtree.r"); cat("parse OK\n")'
```

## 完整运行命令

建议将日志写入 coldmouse 项目的 `logs/` 目录：

```bash
mkdir -p /mnt/disk1/qiuzerui/expriments/coldmouse/logs

/home/zerui/miniconda3/envs/cmintro/bin/Rscript \
  "/home/zerui/code/coldmouse/coldmouse_analysis/using/pseudotime/Pseudotime Trajectory Analysis_monocle2_ddrtree.r" \
  > /mnt/disk1/qiuzerui/expriments/coldmouse/logs/monocle2_ddrtree_pseudotime.log 2>&1
```

## 运行后检查

查看日志末尾：

```bash
tail -80 /mnt/disk1/qiuzerui/expriments/coldmouse/logs/monocle2_ddrtree_pseudotime.log
```

确认输出：

```bash
ls -lh /mnt/disk1/qiuzerui/expriments/coldmouse/files/script_06_monocle2_ddrtree_*
ls -lh /mnt/disk1/qiuzerui/expriments/coldmouse/pictures/monocle2_ddrtree_*
```

成功结束时日志应包含：

```text
Monocle2 DDRTree analysis completed successfully!
```

## 注意事项

1. DDRTree 图是 monocle2 的 component 坐标，不等同于 Seurat UMAP。
2. 如果需要在 UMAP 上展示 pseudotime，请使用 `Pseudotime Trajectory Analysis_monocle2_umap.r`。
3. 如果出现 `dplyr` 或 `igraph` 兼容问题，优先检查脚本开头的 monocle2 compatibility patch。
4. 不要把 `reduceDimension(method = "DDRTree")` 改成 UMAP。monocle2 的 trajectory graph 依赖 DDRTree。
