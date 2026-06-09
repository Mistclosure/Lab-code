# Pseudo-cell CiliaScore

This folder implements a pseudo-cell-level cilia activity scoring method adapted from background-corrected module scoring strategies.

The method is designed for scRNA-seq raw counts. It first aggregates cells into pseudo-cells or pseudo-bulk profiles, then calculates:

```text
CiliaScore = target_score - background_score
```

`target_score` is the average log-normalized expression of cilia genes. `background_score` is the average expression of random control genes matched to cilia genes by both pseudo-cell mean expression and pseudo-cell detection rate.

## Files

```text
R/pseudo_cilia_score.R
scripts/run_pseudo_cilia_score_example.R
README_pseudo_cilia_score.md
```

## Input

`cal_pseudoCiliaScore()` accepts:

- A Seurat object, matrix, or `Matrix::dgCMatrix`.
- Rows must be genes and columns must be cells for matrix input.
- If input is a Seurat object, counts are extracted from `assay` and `layer`.
- If input is a matrix or sparse matrix, `meta` must be supplied.
- Metadata rownames must match cell names. If not identical, overlapping cells are used with a warning.

Required metadata columns:

- `sample_col`: sample or patient ID.
- `group_col`: biological group, such as Primary or Metastasis.

Optional metadata columns:

- `celltype_col`: cell type annotation.
- `subcluster_col`: cell state or malignant subcluster.

## Main Usage

```r
source("R/pseudo_cilia_score.R")

result <- cal_pseudoCiliaScore(
  object = seurat_obj,
  cilia_genes = cilia_genes,
  assay = "RNA",
  layer = "counts",
  sample_col = "orig.ident",
  group_col = "metastasis",
  celltype_col = NULL,
  subcluster_col = "malignant_subcluster",
  target_celltype = NULL,
  pseudo_mode = "chunk",
  cells_per_pseudo = 80,
  min_cells_per_pseudo = 30,
  num.rounds = 1000,
  seed = 123,
  return_matrices = FALSE
)

save_pseudoCiliaScore_results(result, "results/pseudo_cilia_score")
```

## Key Parameters

- `pseudo_mode = "chunk"`: random chunks within each sample/group/cell type/subcluster stratum.
- `pseudo_mode = "sample"`: one pseudo-cell per stratum.
- `cells_per_pseudo = 80`: target cells per chunk.
- `min_cells_per_pseudo = 30`: strata or remainders below this are dropped.
- `min_total_count = 10`: minimum total pseudo-cell counts per gene.
- `min_detect_pseudocells = 3`: minimum number of pseudo-cells with count > 0.
- `min_detect_rate = 0.01`: minimum pseudo-cell detection rate.
- `remove_mt = TRUE`: remove mitochondrial genes matching `^MT-` or `^mt-`.
- `remove_ribo = FALSE`: optionally remove ribosomal genes matching `^RPL|^RPS|^Rpl|^Rps`.
- `n.mean.bins = 20`: bins for mean expression matching.
- `n.detect.bins = 10`: bins for detection-rate matching.
- `n.control = 50`: control genes sampled for each cilia gene in each round.
- `num.rounds = 1000`: random background sampling rounds.
- `case_insensitive_genes = FALSE`: default keeps exact gene-symbol matching.

## Output

The returned object is a list containing:

- `score`: pseudo-cell-level CiliaScore, target score, background score, sample/group labels, and `n_cells`.
- `sample_summary`: sample-level mean, median, SD, pseudo-cell count, and total cell count.
- `pseudo_meta`: metadata for pseudo-cells.
- `pseudo_counts`: filtered pseudo-cell counts if `return_matrices = TRUE`.
- `pseudo_norm`: log-normalized pseudo-cell expression if `return_matrices = TRUE`.
- `gene_info`: mean expression, detection rate, bins, and cilia gene flags.
- `cilia_genes_used`: cilia genes used after matching and filtering.
- `cilia_genes_missing`: missing or filtered-out cilia genes.
- `background_mapping`: matched background pools for each cilia gene.
- `parameters`: analysis parameters.

`save_pseudoCiliaScore_results()` writes:

```text
score.csv
sample_summary.csv
gene_info.csv
cilia_genes_used.txt
cilia_genes_missing.txt
background_mapping.rds
result.rds
```

## Recommended Statistics

Pseudo-cell scores can be plotted to show within-sample or within-state variability, but pseudo-cell chunks should not be interpreted as fully independent patients.

For formal Primary vs Metastasis comparison, use `sample_summary`, where each sample contributes one mean CiliaScore per analyzed stratum.

Example:

```r
wilcox.test(mean_CiliaScore ~ group, data = result$sample_summary)
```

If multiple pseudo-cell chunks are modeled directly, use a mixed model with sample as a random effect:

```r
# Optional example; requires lme4.
# lme4::lmer(CiliaScore ~ group + (1 | sample_id), data = result$score)
```

## Methodological Positioning

This should be described as a pseudo-cell-level cilia activity scoring method adapted from background-corrected module scoring strategies.

Key points:

- It does not score individual cells directly.
- It aggregates cells within sample/cell-state strata to reduce dropout and technical noise.
- Background genes are matched by both mean expression and detection rate.
- Multiple random control rounds improve score stability.
- Higher `target_score - background_score` indicates higher inferred cilia activity.
