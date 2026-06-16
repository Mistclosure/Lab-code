# CRC 单细胞 CRDscore / CiliaScore 随机阴性对照分析执行规划

## 1. 本次改造目标

在 `/home/zerui/code/CRC/cillascore` 下整理并实现新的单细胞层面 CRDscore / CiliaScore 分析流程。该流程从现有五个单数据集 CRDscore 主脚本中保留 CRDscore 计算与绘图相关逻辑，删除 AddModuleScore、FindMarkers、DotPlot、VlnPlot、保存修改后 Seurat 对象等旁支内容，并新增随机同数量基因集阴性对照分析。

明确忽略合并数据集脚本：

- `/home/zerui/code/CRC/GSE132465_GSE231559/03.CRDscore.r`
- `/home/zerui/code/CRC/cholesterol_CRDscore_scripts/GSE132465_GSE231559/03.CRDscore_cholesterol.r`

需要覆盖的五个独立数据集：

- GSE132465
- GSE178318
- GSE178341
- GSE225857
- GSE231559

最终推荐生成两个代码文件：

- `/home/zerui/code/CRC/cillascore/R/singlecell_crdscore_random_control.R`
- `/home/zerui/code/CRC/cillascore/scripts/run_singlecell_crdscore_random_control.R`

结果输出到：

- `/home/zerui/code/CRC/cillascore/results/singlecell_crdscore_random_control/`

该输出目录是本分析专用目录，可以创建新的 `files/` 和 `plots/` 子目录；不要写回或覆盖原始数据集目录中的已有结果。

## 2. 当前项目中可参考的五个主脚本

优先参考单数据集 CRDscore 主脚本，而不是 combined 脚本。

### 2.1 GSE132465

现有脚本：

- `/home/zerui/code/CRC/GSE132465/R/CRDscore_GSE132465.r`

关键逻辑：

- 工作目录：`/mnt/disk1/qiuzerui/downloads/CRC/GSE132465`
- Seurat 对象：`qs/Seurat/Malignant_RNA_assay.qs`
- 临床信息：`metadata/GSE132465_Cli.csv`
- metadata 与临床合并：`meta$orig.ident` 对 `cli$Tumor`
- 分组构造：根据 `TNM stage`，`M0$` 作为 `Primary`，其他作为 `Metastasis`
- 当前基因集示例：`/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/table2_primary_cilia_single_gene_confirmed.csv`
- 当前打分：`CRDscore::cal_CRDscore(expr = data_log2, n.bins = 50, circadians = target_genes, study.type = "scRNAseq")`
- 当前绘图：Stage boxplot 和 metastasis boxplot

改造时只保留：

- 读取对象
- NormalizeData
- 提取 `RNA` assay 的 `data` layer
- 转换到 log2 尺度
- 基因集匹配
- CiliaScore/CRDscore 计算
- Primary vs Metastasis 统计和绘图

删除或避免：

- Stage 分期图不是本次核心，可不迁移
- 不保存修改后的 Seurat 对象
- 不写回原始 GSE132465 目录

### 2.2 GSE178318

现有脚本：

- `/home/zerui/code/CRC/CRDscore_GSE178318.r`

关键逻辑：

- 工作目录：`/mnt/disk1/qiuzerui/downloads/CRC/GSE178318/`
- Seurat 对象：`Malignant.qs`
- 临床信息：`GSE178318_Cli.csv`
- metadata 与临床合并：`meta$orig.ident` 对 `cli$Patient ID`
- 分组构造：现有脚本用 `pbmc1$type <- sub(".*_", "", colnames(pbmc1))`，并保留 `CRC` 和 `LM`
- Primary/Metastasis 标准化建议：`CRC -> Primary`，`LM -> Metastasis`
- 当前代码中 AddModuleScore 已注释，但后续 AddModuleScore plot 仍残留，应删除

注意：

- 现有代码中 `score <- cal_CRDscore(expr = exp, ...)`，其中 `exp = data_log2`
- 需要确认 `colnames(pbmc1)` 中最后一个 `_` 后的字段是否稳定为 `CRC`/`LM`

### 2.3 GSE178341

现有脚本：

- `/home/zerui/code/CRC/CRDscore_GSE178341.r`

关键逻辑：

- 工作目录：`/mnt/disk1/qiuzerui/downloads/CRC/GSE178341/`
- Seurat 对象：`Malignant.qs`
- 临床信息：`Cli.csv`
- metadata 与临床合并：`meta$orig.ident` 对 `cli$PatientBarcode`
- 现有分组字段：`Tumor Stage`

风险点：

- 用户目标要求至少能区分 `Primary` 和 `Metastasis`。现有脚本使用 `Tumor Stage`，不一定是转移状态。
- 执行时必须在 metadata 和 clinical 中自动搜索候选字段，例如 `metastasis`、`Metastasis`、`type`、`Tumor Stage`、`Stage`，但只有能明确映射到 Primary/Metastasis 时才运行。
- 如果 GSE178341 无法可靠区分 Primary/Metastasis，应 warning 并跳过该数据集，不中断其他数据集。

### 2.4 GSE225857

现有脚本：

- `/home/zerui/code/CRC/CRDscore_GSE225857.r`

关键逻辑：

- 工作目录：`/mnt/disk1/qiuzerui/downloads/CRC/GSE225857/`
- Seurat 对象：`Malignant.qs`
- 临床信息：`GSE225857_Cli.csv`
- metadata 与临床合并：`meta$orig.ident` 对 `cli$Patient ID`
- 分组字段：`Liver metastasis (LM)`
- 分组构造：`Yes -> Metastasis`，其他或 `No -> Primary`

删除内容：

- AddModuleScore
- VlnPlot
- FindMarkers
- DotPlot
- top genes 相关代码

### 2.5 GSE231559

现有脚本：

- `/home/zerui/code/CRC/GSE231559/CRDscore_GSE231559.r`

关键逻辑：

- 工作目录：`/mnt/disk1/qiuzerui/downloads/CRC/GSE231559`
- Seurat 对象：`Malignant_RNA_assay.qs`
- 临床信息：`GSE231559_Cli.csv`
- metadata 与临床合并：现有脚本使用 `merge(meta, cli, by = "Patient")`
- 分组字段候选：`metastasis`
- 当前脚本会执行：
  - `pbmc1[['RNA']]$data_log2 = data_log2`
  - `qsave(pbmc1, 'Malignant_GSE231559.qs')`

必须删除或避免：

- 不要修改 Seurat object
- 不要保存 `Malignant_GSE231559.qs`
- 只在内存中使用 `data_log2`

## 3. 推荐的新项目结构

在 `/home/zerui/code/CRC/cillascore` 下创建：

```text
cillascore/
  R/
    singlecell_crdscore_random_control.R
  scripts/
    run_singlecell_crdscore_random_control.R
  results/
    singlecell_crdscore_random_control/
      files/
      plots/
      logs/
```

其中：

- `R/singlecell_crdscore_random_control.R`：放通用函数。
- `scripts/run_singlecell_crdscore_random_control.R`：放配置、数据集循环、结果汇总、绘图调用。
- `results/.../files/`：保存 CSV。
- `results/.../plots/`：保存 PDF/PNG。
- `results/.../logs/`：保存 warning、跳过原因、sessionInfo。

## 4. 输入配置设计

在主脚本中定义一个 `dataset_configs` list，每个数据集一个配置，避免五个脚本复制粘贴。

建议字段：

```r
dataset_configs <- list(
  GSE132465 = list(
    dataset = "GSE132465",
    work_dir = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465",
    object_path = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/qs/Seurat/Malignant_RNA_assay.qs",
    clinical_path = "/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/metadata/GSE132465_Cli.csv",
    object_reader = "qread",
    assay = "RNA",
    layer = "data",
    sample_col_meta = "orig.ident",
    clinical_merge_by_meta = "orig.ident",
    clinical_merge_by_clinical = "Tumor",
    group_source = "TNM stage",
    group_rule = "tnm_m0_primary"
  )
)
```

其余数据集配置：

```text
GSE178318:
  object_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE178318/Malignant.qs
  clinical_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE178318/GSE178318_Cli.csv
  merge: orig.ident -> Patient ID
  group_rule = cellname_suffix_crc_lm

GSE178341:
  object_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE178341/Malignant.qs
  clinical_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE178341/Cli.csv
  merge: orig.ident -> PatientBarcode
  group_rule = auto_detect_primary_metastasis
  note: 若无法明确 Primary/Metastasis，则 warning 并跳过

GSE225857:
  object_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE225857/Malignant.qs
  clinical_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE225857/GSE225857_Cli.csv
  merge: orig.ident -> Patient ID
  group_source = Liver metastasis (LM)
  group_rule = yes_no_liver_metastasis

GSE231559:
  object_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE231559/Malignant_RNA_assay.qs
  clinical_path = /mnt/disk1/qiuzerui/downloads/CRC/GSE231559/GSE231559_Cli.csv
  merge: Patient -> Patient
  group_source = metastasis
  group_rule = metastasis_column
```

所有路径在运行前用 `file.exists()` 检查。若对象或临床文件不存在：

```r
warning(sprintf("[%s] input path not found, skip dataset: %s", dataset, path))
return(NULL)
```

不要因为单个数据集缺失而停止全流程。

## 5. cilia gene set 输入设计

用户要求使用人工整理的 cilia gene set xlsx，包含 `Sheet1` 和 `Sheet2`。

主脚本应定义：

```r
gene_set_xlsx <- "/path/to/cilia_gene_list.xlsx"
gene_set_sheets <- c("Sheet1", "Sheet2")
gene_column_candidates <- c("Gene", "gene", "genes", "symbol", "Symbol", "gene_symbol", "GeneSymbol")
```

如果当前项目中没有指定 xlsx，可先按以下顺序自动查找：

1. 用户显式传入的 `gene_set_xlsx`
2. `/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/input/` 下含 `cilia`、`primary_cilia`、`single_gene` 的 xlsx
3. `/home/zerui/code/CRC/results/cilia_gene_score_correlation/all_gene_consistency_summary.xlsx`
4. 如果仍找不到，停止并给出清晰报错，要求用户提供 xlsx 路径

注意：

- 本规划不建议继续使用现有 CSV 作为唯一输入，因为用户明确要求 xlsx 且含 Sheet1/Sheet2。
- 可兼容 CSV 作为 fallback，但正式流程应优先 xlsx。

读取函数 `load_cilia_gene_sets()` 应返回：

```r
list(
  Sheet1 = c("GENE1", "GENE2", ...),
  Sheet2 = c("GENE3", "GENE4", ...)
)
```

并保存每个数据集、每个 gene set 的 missing genes：

```text
files/{dataset}_{gene_set_name}_missing_genes.csv
```

字段：

- dataset
- gene_set_name
- gene
- status: missing_from_expression

若真实 gene set 与表达矩阵交集少于 5 个基因，跳过该 gene set 并 warning。

## 6. 表达矩阵处理要求

函数 `extract_log_normalized_matrix(seurat_obj, assay = "RNA", layer = "data")`：

1. 不修改原始 Seurat 对象。
2. 若 `RNA/data` layer 不存在，则在对象副本或临时变量中调用 `NormalizeData()`。
3. 提取 log-normalized expression。
4. 确保行为基因、列为细胞。
5. 过滤全 0 或不可检测基因：
   - 若 counts layer 存在，用 `Matrix::rowSums(counts) > 0`
   - 若 counts 不存在，用 `Matrix::rowSums(data > 0) > 0`
6. 对 CRDscore 包兼容时，现有脚本使用 `LayerData(..., layer = "data") / log(2)` 将自然 log 转为 log2。新流程可保留这一点，保证与现有 CRDscore 结果方向和尺度接近。

建议不要把大型稀疏矩阵强制转成 `data.frame`，除非 `CRDscore::cal_CRDscore()` 必须要求。优先保持 `matrix` 或 `dgCMatrix`，降低内存压力。

## 7. 单细胞 CiliaScore / CRDscore 计算方法

优先复用现有 `CRDscore::cal_CRDscore()`，但它不一定返回 target/background 分项。用户要求输出：

- `CiliaScore`
- `target_score`
- `background_score`

因此建议实现一个本项目内部的背景校正单细胞打分函数，命名为：

```r
calculate_background_corrected_score()
```

算法：

1. 输入 log-normalized expression 矩阵 `expr`，行为基因、列为细胞。
2. 计算每个基因平均表达：
   ```r
   gene_avg <- Matrix::rowMeans(expr)
   ```
3. 按平均表达分 bin：
   ```r
   bins <- cut(rank(gene_avg, ties.method = "random"), breaks = n_bins, labels = FALSE)
   ```
4. 对目标基因集取表达矩阵中存在的基因。
5. 对每个目标基因，在同一个表达 bin 中抽取背景基因：
   - 排除真实目标基因
   - 每个目标基因抽 `ctrl_per_gene` 个，默认 25 或 50
   - 如果同 bin 候选不足，可允许放回抽样，或扩展到相邻 bin
6. 每个细胞计算：
   ```r
   target_score <- Matrix::colMeans(expr[target_genes_used, , drop = FALSE])
   background_score <- Matrix::colMeans(expr[background_genes, , drop = FALSE])
   CiliaScore <- target_score - background_score
   ```
7. 返回：
   ```r
   list(
     score_df = data.frame(
       cell_id,
       CiliaScore,
       target_score,
       background_score
     ),
     target_genes_used = target_genes_used,
     background_genes_used = background_genes_used,
     bins = bins
   )
   ```

说明：

- 该实现符合用户要求的 `target_gene_score - background_gene_score`。
- 如果需要和旧脚本结果对齐，可额外保存 `CRDscore_package_score`，但不要作为主输出字段，避免结果表结构变复杂。
- 旧脚本中的 `cal_CRDscore()` 可作为备选检查，不作为唯一实现。

## 8. 随机阴性对照设计

主参数：

```r
n_random <- 1000
seed <- 12345
random_modes <- c("simple_random", "expression_matched_random")
expected_direction <- "Metastasis_high"
```

支持扩展：

```r
n_random <- 10000
```

但默认 1000，避免五个数据集、两个 sheet、两种随机模式运行过慢。

### 8.1 候选基因池

候选基因为表达矩阵中的可检测基因：

```r
detected_genes <- rownames(expr)[Matrix::rowSums(expr > 0) > 0]
candidate_genes <- setdiff(detected_genes, all_real_cilia_genes_present_or_all_input_cilia_genes)
```

注意：

- 随机基因不能包含真实 cilia genes。
- 推荐从当前数据集可检测基因中抽取。
- 对 Sheet1 和 Sheet2 分别运行时，排除对应 sheet 的真实基因；更严格做法是排除 Sheet1 + Sheet2 的并集，推荐采用严格做法。

### 8.2 simple_random

每次随机：

```r
random_genes <- sample(candidate_genes, size = n_target_genes, replace = FALSE)
```

如果候选基因不足：

- warning
- 该 gene set / random mode 跳过

### 8.3 expression_matched_random

每次随机：

1. 对真实 target genes 获取表达 bin。
2. 对每个真实 target gene，从同 bin 的候选基因中抽 1 个。
3. 排除真实 cilia genes。
4. 尽量不重复。
5. 如果同 bin 候选不足：
   - 先尝试相邻 bin，例如 `bin - 1`、`bin + 1`
   - 仍不足时允许放回抽样
   - 记录 warning 但不停止

该模式的随机 gene set 与真实 gene set 平均表达分布更接近，是主要阴性对照。

### 8.4 随机 gene set 不存在基因处理

因为随机候选来自表达矩阵 rownames，正常不会缺失。若后续逻辑中发现可用基因数小于真实 gene set 的 80% 或少于 5：

- 重新抽样，最多尝试 `max_resample_attempts = 20`
- 仍失败则跳过该 random_id 并 warning

## 9. Primary vs Metastasis 统计定义

统一分组为：

```text
Primary
Metastasis
```

所有数据集输出前必须映射到这两个标准标签。

### 9.1 delta mean

定义：

```r
delta_mean <- mean(score[group == "Metastasis"]) - mean(score[group == "Primary"])
```

即正值表示 Metastasis 更高。

### 9.2 delta median

```r
delta_median <- median(score[group == "Metastasis"]) - median(score[group == "Primary"])
```

### 9.3 p value 和 neg_log10_p

使用 Wilcoxon rank-sum test：

```r
p_value <- wilcox.test(score ~ group)$p.value
neg_log10_p <- -log10(pmax(p_value, .Machine$double.xmin))
```

不要比较形如 `< 2.2e-16` 的打印字符串，必须比较数值化后的 `neg_log10_p`。

### 9.4 effect size

建议同时保留一个明确命名的效应量：

```r
effect_size <- delta_mean / pooled_sd
```

其中：

```r
pooled_sd <- sqrt(((n_primary - 1) * sd_primary^2 + (n_metastasis - 1) * sd_metastasis^2) / (n_primary + n_metastasis - 2))
```

若 `pooled_sd == 0`，返回 `NA_real_`。

### 9.5 AUC

定义为随机抽取一个 Metastasis 细胞得分大于一个 Primary 细胞得分的概率。

可使用 `pROC::auc()` 或自己用 Wilcoxon U 统计量计算：

```r
rank_sum_meta <- sum(rank(score)[group == "Metastasis"])
u_meta <- rank_sum_meta - n_meta * (n_meta + 1) / 2
auc <- u_meta / (n_meta * n_primary)
```

解释：

- `AUC > 0.5` 表示 Metastasis 倾向更高。
- `AUC < 0.5` 表示 Primary 倾向更高。

### 9.6 Cliff's delta

用户写的是 “AUC 或 Cliff's delta”。建议主输出字段保留 `AUC`，因为输出模板指定了 AUC。若想补充，Cliff's delta 可由 AUC 转换：

```r
cliffs_delta <- 2 * auc - 1
```

但为了不改变用户指定字段，`random_null_results.csv` 和 `observed_group_test.csv` 中只必需写 `AUC` 和 `effect_size`。

### 9.7 direction

根据 `delta_mean`：

```r
direction <- dplyr::case_when(
  delta_mean > 0 ~ "Metastasis_high",
  delta_mean < 0 ~ "Primary_high",
  TRUE ~ "No_difference"
)
```

## 10. 经验 P 值和 percentile 计算

对每个 dataset、gene_set_name、random_mode 分别计算。

默认预期方向：

```r
expected_direction <- "Metastasis_high"
```

单侧经验 P：

```r
empirical_p_by_delta <- (sum(random_delta_mean >= observed_delta_mean) + 1) / (n_random + 1)
empirical_p_by_neg_log10_p <- (sum(random_neg_log10_p >= observed_neg_log10_p) + 1) / (n_random + 1)
empirical_p_by_effect_size <- (sum(random_effect_size >= observed_effect_size) + 1) / (n_random + 1)
```

若不预设方向，则用双侧：

```r
empirical_p_abs_delta <- (sum(abs(random_delta_mean) >= abs(observed_delta_mean)) + 1) / (n_random + 1)
empirical_p_abs_effect_size <- (sum(abs(random_effect_size) >= abs(observed_effect_size)) + 1) / (n_random + 1)
```

percentile 推荐：

```r
percentile_by_delta <- mean(random_delta_mean <= observed_delta_mean, na.rm = TRUE)
percentile_by_neg_log10_p <- mean(random_neg_log10_p <= observed_neg_log10_p, na.rm = TRUE)
percentile_by_effect_size <- mean(random_effect_size <= observed_effect_size, na.rm = TRUE)
```

如果最终采用双侧解释，应额外说明 percentile 仍是按原始方向给出的，经验 P 使用绝对值。

## 11. 输出文件设计

每个 dataset / gene_set_name / random_mode 保存独立文件，同时保存跨数据集汇总文件。

### 11.1 observed_score.csv

建议文件名：

```text
files/{dataset}_{gene_set_name}_observed_score.csv
```

字段：

- cell_id
- sample_id
- group
- CiliaScore
- target_score
- background_score

`sample_id` 从标准化后的 metadata 字段中提取，优先：

1. `orig.ident`
2. `Patient`
3. `sample`
4. 如果不存在，则填 `NA`

### 11.2 observed_group_test.csv

建议文件名：

```text
files/{dataset}_{gene_set_name}_observed_group_test.csv
```

字段：

- dataset
- gene_set_name
- n_genes_used
- mean_primary
- mean_metastasis
- delta_mean
- median_primary
- median_metastasis
- delta_median
- p_value
- neg_log10_p
- effect_size
- AUC
- direction

同时保存跨数据集汇总：

```text
files/all_observed_group_test.csv
```

### 11.3 random_null_results.csv

建议文件名：

```text
files/{dataset}_{gene_set_name}_{random_mode}_random_null_results.csv
```

字段：

- dataset
- gene_set_name
- random_mode
- random_id
- n_genes
- delta_mean
- p_value
- neg_log10_p
- effect_size
- AUC
- direction

同时保存跨数据集汇总：

```text
files/all_random_null_results.csv
```

### 11.4 empirical_test_summary.csv

建议文件名：

```text
files/{dataset}_{gene_set_name}_{random_mode}_empirical_test_summary.csv
```

字段：

- dataset
- gene_set_name
- random_mode
- observed_delta_mean
- observed_p_value
- observed_neg_log10_p
- observed_effect_size
- observed_AUC
- random_delta_mean_mean
- random_delta_mean_sd
- random_neg_log10_p_mean
- random_neg_log10_p_sd
- percentile_by_delta
- percentile_by_neg_log10_p
- percentile_by_effect_size
- empirical_p_by_delta
- empirical_p_by_neg_log10_p
- empirical_p_by_effect_size

同时保存跨数据集汇总：

```text
files/all_empirical_test_summary.csv
```

### 11.5 日志与跳过记录

保存：

```text
logs/skipped_datasets_or_gene_sets.csv
logs/sessionInfo.txt
```

`skipped_datasets_or_gene_sets.csv` 字段：

- dataset
- gene_set_name
- step
- reason

常见 reason：

- object_path_not_found
- clinical_path_not_found
- no_primary_metastasis_group
- fewer_than_5_genes_overlap
- random_candidate_pool_too_small

## 12. 绘图设计

所有图保存 PDF 和 PNG。绘图用 `ggplot2`，publication style。不要使用旧脚本中的固定 `ylim(-0.2, 0.2)`，避免截断不同数据集的真实分布。

### 12.1 随机 delta mean null distribution

每个 dataset / gene_set / random_mode 一张图：

```text
plots/{dataset}_{gene_set_name}_{random_mode}_null_delta_mean.pdf
plots/{dataset}_{gene_set_name}_{random_mode}_null_delta_mean.png
```

图形：

- x = `random_delta_mean`
- histogram 或 density
- 竖线标注 `observed_delta_mean`
- 标题包含 dataset、gene_set_name、random_mode、empirical P

### 12.2 随机 -log10(p) distribution

文件：

```text
plots/{dataset}_{gene_set_name}_{random_mode}_null_neg_log10_p.pdf
plots/{dataset}_{gene_set_name}_{random_mode}_null_neg_log10_p.png
```

图形：

- x = `random_neg_log10_p`
- 竖线标注 `observed_neg_log10_p`

### 12.3 observed vs random summary plot

跨数据集汇总图：

```text
plots/observed_vs_random_delta_mean_summary.pdf
plots/observed_vs_random_delta_mean_summary.png
plots/observed_vs_random_effect_size_summary.pdf
plots/observed_vs_random_effect_size_summary.png
```

图形：

- y = dataset
- x = delta_mean 或 effect_size
- 灰色 violin / boxplot / jitter 表示 random
- 彩色点表示 observed
- facet by gene_set_name 和 random_mode

### 12.4 每个 dataset 的 CiliaScore boxplot

文件：

```text
plots/{dataset}_{gene_set_name}_observed_CiliaScore_boxplot.pdf
plots/{dataset}_{gene_set_name}_observed_CiliaScore_boxplot.png
```

图形：

- x = `group`
- y = `CiliaScore`
- group 顺序固定为 `Primary`, `Metastasis`
- 显示 Wilcoxon p value、delta_mean、AUC
- 点太多时不画所有散点，最多使用 sample 或只画 boxplot + violin，避免文件过大

## 13. 通用函数拆分建议

放入 `/home/zerui/code/CRC/cillascore/R/singlecell_crdscore_random_control.R`。

建议函数：

```r
safe_dir_create(path)
load_required_packages()
load_cilia_gene_sets(gene_set_xlsx, sheets, gene_column_candidates)
load_seurat_object(path, reader = c("qread", "readRDS"))
extract_log_normalized_matrix(seurat_obj, assay = "RNA", layer = "data", log2_transform = TRUE)
standardize_dataset_metadata(seurat_obj, clinical_df, config)
derive_primary_metastasis_group(meta_clinical_df, config)
match_gene_set(genes, expr_genes, dataset, gene_set_name)
make_expression_bins(expr, n_bins = 50, seed = 12345)
sample_background_genes(target_genes, gene_bins, candidate_genes, ctrl_per_gene = 25, seed = NULL)
calculate_background_corrected_score(expr, target_genes, gene_bins, candidate_genes, ctrl_per_gene = 25, seed = NULL)
sample_random_gene_set_simple(candidate_genes, n_genes, seed = NULL)
sample_random_gene_set_expression_matched(target_genes, gene_bins, candidate_genes, seed = NULL)
calculate_group_statistics(score_df, dataset, gene_set_name)
calculate_auc(score, group)
run_random_null_for_gene_set(expr, observed_genes, all_real_genes, group_df, gene_bins, n_random, random_mode, seed)
summarize_empirical_test(observed_stats, random_stats, expected_direction = "Metastasis_high")
plot_observed_boxplot(observed_score_df, observed_stats, output_prefix)
plot_random_delta_distribution(random_df, observed_stats, empirical_summary, output_prefix)
plot_random_neglogp_distribution(random_df, observed_stats, empirical_summary, output_prefix)
plot_observed_vs_random_summary(all_observed, all_random, output_dir)
write_skip_log(skip_records, output_path)
```

函数返回值尽量使用 list 或 tibble，避免依赖全局变量。

## 14. 主脚本执行流程

`scripts/run_singlecell_crdscore_random_control.R` 结构建议：

```r
PROJECT_DIR <- "/home/zerui/code/CRC/cillascore"
source(file.path(PROJECT_DIR, "R", "singlecell_crdscore_random_control.R"))

load_required_packages()

set.seed(12345)

params <- list(
  n_random = 1000,
  n_bins = 50,
  ctrl_per_gene = 25,
  seed = 12345,
  random_modes = c("simple_random", "expression_matched_random"),
  expected_direction = "Metastasis_high"
)

output_dir <- file.path(PROJECT_DIR, "results", "singlecell_crdscore_random_control")
files_dir <- file.path(output_dir, "files")
plots_dir <- file.path(output_dir, "plots")
logs_dir <- file.path(output_dir, "logs")

gene_sets <- load_cilia_gene_sets(...)

all_observed_scores <- list()
all_observed_tests <- list()
all_random_results <- list()
all_empirical_summaries <- list()
skip_records <- list()

for (dataset_name in names(dataset_configs)) {
  result <- tryCatch(
    run_one_dataset(...),
    error = function(e) {
      warning(sprintf("[%s] failed: %s", dataset_name, conditionMessage(e)))
      NULL
    }
  )
}

write_csv(bind_rows(all_observed_tests), file.path(files_dir, "all_observed_group_test.csv"))
write_csv(bind_rows(all_random_results), file.path(files_dir, "all_random_null_results.csv"))
write_csv(bind_rows(all_empirical_summaries), file.path(files_dir, "all_empirical_test_summary.csv"))

plot_observed_vs_random_summary(...)
writeLines(capture.output(sessionInfo()), file.path(logs_dir, "sessionInfo.txt"))
```

核心要求：

- 单个数据集失败不影响其他数据集。
- 单个 gene set 失败不影响同一数据集其他 gene set。
- 单个 random mode 失败不影响 observed score 输出。
- 任何跳过都要写入 skip log。

## 15. 需要从旧脚本删除的内容

如果后续要将五个脚本也复制整理到 `/home/zerui/code/CRC/cillascore`，应删除：

- `AddModuleScore(...)`
- `AddModuleScore1` 相关 metadata 提取和绘图
- `VlnPlot(...)`
- `FindMarkers(...)`
- `DotPlot(...)`
- top 30 gene plotting
- Stage 分期图，除非用户额外要求
- `pbmc1[['RNA']]$data_log2 = data_log2`
- `qsave(pbmc1, ...)`
- 写入原始数据集目录的 CSV/PNG
- 全局 `setwd()` 依赖；改为显式绝对路径

保留：

- 各数据集对象路径
- 临床文件路径
- metadata 与临床表合并规则
- Primary/Metastasis 构造规则
- `NormalizeData()` 与 `LayerData(..., layer = "data")`
- CRDscore/CiliaScore boxplot 逻辑

## 16. 性能与内存优化

五个单细胞数据集加 1000 次随机可能较慢。建议：

1. 不把整个表达矩阵转成 dense `data.frame`，尽量保持 `dgCMatrix`。
2. 预先计算并缓存：
   - `gene_avg`
   - `gene_bins`
   - `candidate_genes`
   - `group` 向量
3. 随机 score 只保留统计结果，不保存每次随机的逐细胞 score。
4. observed score 才保存逐细胞 CSV。
5. 每完成一个 dataset 后执行 `gc()`。
6. 可选支持 `future.apply` 或 `BiocParallel` 并行，但第一版建议串行，保证结果稳定和 debug 简单。
7. 若 `n_random = 10000`，建议只开启 `expression_matched_random` 或按数据集分批运行。

## 17. 鲁棒性检查清单

实现时必须检查：

- 输入对象路径是否存在。
- 临床文件路径是否存在。
- Seurat 对象是否可读取。
- `RNA` assay 是否存在。
- `data` layer 是否存在；不存在则 NormalizeData。
- 表达矩阵是否有行名和列名。
- metadata 行名是否对应细胞 ID。
- clinical merge 后是否保留足够细胞。
- group 是否同时包含 `Primary` 和 `Metastasis`。
- 每组细胞数是否大于 1。
- 每个 gene set 与表达矩阵交集是否至少 5 个基因。
- 随机候选基因池是否足够。
- `p_value` 是否为数值且非 NA。
- `neg_log10_p` 使用 `pmax(p_value, .Machine$double.xmin)`。
- 输出文件名不含空格或特殊字符，gene_set_name 用安全版本，例如 `make.names()` 或自定义 `sanitize_name()`。

## 18. 最终结果解释模板

运行完成后，应基于 `all_observed_group_test.csv` 和 `all_empirical_test_summary.csv` 回答用户关心的五个问题。

### 18.1 真实 cilia gene set 的 Primary vs Metastasis 差异是否显著？

判断依据：

- `observed_p_value < 0.05`
- `observed_neg_log10_p`
- `delta_mean` 和 `AUC` 的方向

注意说明：

- 单细胞细胞数大时 p value 容易极小。
- 不能只看 p value，必须同时看 effect size、delta_mean、AUC。

### 18.2 真实 gene set 的显著性是否超过随机同数量 gene set？

判断依据：

- `percentile_by_neg_log10_p` 是否接近 1，例如 > 0.95
- `empirical_p_by_neg_log10_p < 0.05`

### 18.3 真实 gene set 的 effect size 是否超过随机同数量 gene set？

判断依据：

- `percentile_by_effect_size > 0.95`
- `empirical_p_by_effect_size < 0.05`

### 18.4 哪些数据集中真实 cilia gene set 明显优于随机阴性对照？

推荐定义：

一个 dataset / gene_set / random_mode 被认为明显优于随机，需同时满足：

- observed direction 为 `Metastasis_high`，或符合用户指定预期方向
- `empirical_p_by_delta < 0.05`
- `empirical_p_by_neg_log10_p < 0.05`
- `empirical_p_by_effect_size < 0.05`
- `observed_AUC > 0.5`

更严格时可要求 `percentile_by_delta`、`percentile_by_neg_log10_p`、`percentile_by_effect_size` 都 > 0.95。

### 18.5 如果某数据集中真实 cilia gene set 不优于随机

解释模板：

```text
在 {dataset} 的 {gene_set_name} 中，真实 cilia gene set 的 observed delta/effect size 未超过随机同数量 gene set 的主要分布，或 empirical P 未达到显著。这说明该数据集中 cilia signal 在 Primary vs Metastasis 层面不稳定，可能受样本构成、转移标签、肿瘤细胞筛选、表达深度或 gene set overlap 影响。该结论不否定其他数据集中的 cilia signal，需要结合跨数据集一致性判断。
```

## 19. 验收标准

完成后应能在 R 中运行：

```bash
Rscript /home/zerui/code/CRC/cillascore/scripts/run_singlecell_crdscore_random_control.R
```

至少生成：

```text
/home/zerui/code/CRC/cillascore/results/singlecell_crdscore_random_control/files/all_observed_group_test.csv
/home/zerui/code/CRC/cillascore/results/singlecell_crdscore_random_control/files/all_random_null_results.csv
/home/zerui/code/CRC/cillascore/results/singlecell_crdscore_random_control/files/all_empirical_test_summary.csv
/home/zerui/code/CRC/cillascore/results/singlecell_crdscore_random_control/logs/sessionInfo.txt
```

每个成功运行的 dataset / gene_set 至少生成：

- observed_score.csv
- observed_group_test.csv
- random_null_results.csv
- empirical_test_summary.csv
- CiliaScore boxplot
- random delta null distribution plot
- random neg_log10_p null distribution plot

如果某个数据集路径不存在或无法识别 Primary/Metastasis：

- 控制台 warning
- 写入 skip log
- 其他数据集继续运行

## 20. 建议的实现顺序

1. 创建目录结构 `R/`、`scripts/`、`results/.../{files,plots,logs}`。
2. 写通用函数文件，先实现：
   - 包加载
   - gene set 读取
   - Seurat 对象读取
   - 表达矩阵提取
   - metadata/clinical 合并
   - Primary/Metastasis 标准化
3. 先只跑 observed CiliaScore，不做随机。
4. 验证五个数据集中哪些能成功得到 `Primary` 和 `Metastasis`。
5. 加入 `calculate_group_statistics()`，输出 observed_group_test。
6. 加入 `simple_random`，先用 `n_random = 10` 做 smoke test。
7. 加入 `expression_matched_random`，先用 `n_random = 10` 做 smoke test。
8. 确认输出表字段完全符合用户要求。
9. 加入全部绘图。
10. 将 `n_random` 改为 1000 正式运行。
11. 检查 skip log 和 missing genes。
12. 汇总最终结论。

## 21. 重要注意事项

- 不要使用 pseudo-cell。
- 不要使用 bulk 分支。
- 不要合并不同数据集后再打分。
- 每个数据集内部独立计算 score。
- 不要修改原始 Seurat 对象。
- 不要保存修改后的 Seurat 对象。
- 不要继续输出 AddModuleScore 结果。
- 不要只比较 p value；必须比较 `neg_log10_p`、delta_mean、effect size、AUC 和 empirical P。
- GSE178341 目前从已有脚本看只明确使用 `Tumor Stage`，不一定满足 Primary/Metastasis，对该数据集要特别检查，不能强行把 stage 当 metastasis，除非临床表中有明确转移字段或用户确认映射规则。
