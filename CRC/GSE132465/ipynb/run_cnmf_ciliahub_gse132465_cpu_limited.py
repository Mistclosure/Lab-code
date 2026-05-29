# %% [markdown]
# # GSE132465 恶性细胞 CiliaHub 基因 cNMF 分析
# 
# 本 notebook 基于当前项目已经生成的恶性细胞 AnnData 对象继续分析，不从 GSE132465 原始矩阵重新开始。
# 
# 默认输入对象来自已有流程：
# 
# - `/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/Malignant_RNA_assay_for_cNMF_countsX.h5ad`
# - 该对象由 `ipynb/Malignant_to_py.ipynb` / `ipynb/Malignant_to_py.py` 从 `Malignant_RNA_assay.qs` 导出，`X` 和 `layers['counts']` 都保存 raw counts。
# 
# 输出会写入独立目录：
# 
# - `/mnt/disk1/qiuzerui/downloads/CRC/GSE132465/cNMF_CiliaHub_GSE132465/`

# %% [markdown]
# ## 1. 环境和路径设置

# %%
# ==============================================================================
# 1. 环境和路径设置
# ==============================================================================

from pathlib import Path
import os
import re
import shutil
import subprocess
import warnings

# matplotlib / numba 在部分服务器上默认缓存目录不可写，这里必须在导入相关包前固定到 /tmp
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-cnmf")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

# ------------------------------------------------------------------------------
# CPU 线程限制：必须放在 numpy / scipy / sklearn / cnmf 导入之前
# ------------------------------------------------------------------------------
# 你的服务器有 192 cores，但 cNMF/NMF 通常不是核心越多越快；
# 10000 cells x 2000 genes 建议先用 16 核，最多可尝试 24 或 32 核。
# 运行时可用 shell 临时覆盖：CNMF_CPU_THREADS=24 python run_xxx.py
CPU_THREADS = int(os.environ.get("CNMF_CPU_THREADS", "16"))

for _var in [
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "BLIS_NUM_THREADS",
]:
    os.environ.setdefault(_var, str(CPU_THREADS))

# 防止 MKL/OpenMP 动态加线程导致偷偷吃满 CPU
os.environ.setdefault("MKL_DYNAMIC", "FALSE")
os.environ.setdefault("OMP_DYNAMIC", "FALSE")

import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib.pyplot as plt

# 再用 threadpoolctl 做一层运行时限制；如果没安装，不影响主流程
try:
    from threadpoolctl import threadpool_limits, threadpool_info
    _threadpool_controller = threadpool_limits(limits=CPU_THREADS)
    print(f"已限制 BLAS/OpenMP 线程数为: {CPU_THREADS}")
except ImportError:
    threadpool_info = None
    _threadpool_controller = None
    print("提示：未安装 threadpoolctl，仅使用环境变量限制线程。")

# 可视化包不是 cNMF 核心依赖，但后面画 heatmap 会用到
try:
    import seaborn as sns
except ImportError:
    sns = None
    print("警告：未安装 seaborn，后续 heatmap 可视化单元会跳过。")

try:
    import anndata as ad
except ImportError as e:
    raise ImportError(
        "缺少 anndata。请先在当前 Python 环境安装：pip install anndata"
    ) from e

# cNMF 包；如果这里报错，先安装 cnmf 后重启 kernel
try:
    from cnmf import cNMF
except ImportError as e:
    raise ImportError(
        "缺少 cnmf。建议使用 conda 环境安装，见 notebook 末尾说明。"
    ) from e

# ------------------------------------------------------------------------------
# 关键路径：后续如果换输入对象或输出目录，只需要改这里
# ------------------------------------------------------------------------------
WORK_DIR = Path("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465")
INPUT_H5AD = WORK_DIR / "Malignant_RNA_assay_for_cNMF_countsX.h5ad"
CILIAHUB_CSV = WORK_DIR / "ciliahub_genes_list.csv"

# 调试开关：正式分析保持默认 0；需要快速检查全流程时，在 shell 里设置 CNMF_FAST_DEBUG=1。
FAST_DEBUG = os.environ.get("CNMF_FAST_DEBUG", "0") == "1"

OUT_DIR = WORK_DIR / ("cNMF_CiliaHub_GSE132465_debug" if FAST_DEBUG else "cNMF_CiliaHub_GSE132465")
CNMF_NAME = "GSE132465_malignant_CiliaHub_cNMF_debug" if FAST_DEBUG else "GSE132465_malignant_CiliaHub_cNMF"
CNMF_RUN_DIR = OUT_DIR / CNMF_NAME
PLOTS_DIR = OUT_DIR / "plots"
TABLES_DIR = OUT_DIR / "tables"
TMP_DIR = OUT_DIR / "tmp"

for d in [OUT_DIR, PLOTS_DIR, TABLES_DIR, TMP_DIR]:
    d.mkdir(parents=True, exist_ok=True)

# ------------------------------------------------------------------------------
# cNMF 参数：先用相对保守的范围；可根据 K selection plot 再调整
# ------------------------------------------------------------------------------
RANDOM_SEED = 9
if FAST_DEBUG:
    # 快速调试模式只用于检查 notebook 是否能完整跑通；正式分析请不要设置 CNMF_FAST_DEBUG。
    N_ITER = int(os.environ.get("CNMF_DEBUG_N_ITER", "3"))
    K_VALUES = np.arange(2, 4)
else:
    N_ITER = 100
    K_VALUES = np.arange(2, 13)  # CiliaHub 限定基因集通常不需要过大的 K
DENSITY_THRESHOLD = 0.5
LOCAL_NEIGHBORHOOD_SIZE = 0.5 if FAST_DEBUG else 0.3
TOTAL_WORKERS = 1
WORKER_INDEX = 0

# 过滤参数：细胞 QC 先在全基因 raw counts 上做；基因过滤在 CiliaHub 子集上做
MIN_CELL_TOTAL_COUNTS = 500
MIN_CELL_DETECTED_GENES = 100
MIN_GENE_EXPRESSED_CELLS = 20
MIN_GENE_TOTAL_COUNTS = 10

# 如果已经看过 K selection plot，可以手动改这里；None 时默认用 K_VALUES 中间值继续跑 consensus
SELECTED_K = int(K_VALUES[0]) if FAST_DEBUG else None

print("WORK_DIR:", WORK_DIR)
print("INPUT_H5AD:", INPUT_H5AD)
print("CILIAHUB_CSV:", CILIAHUB_CSV)
print("OUT_DIR:", OUT_DIR)
print("FAST_DEBUG:", FAST_DEBUG)
print("K_VALUES:", list(K_VALUES))
print("N_ITER:", N_ITER)
print("DENSITY_THRESHOLD:", DENSITY_THRESHOLD)
print("LOCAL_NEIGHBORHOOD_SIZE:", LOCAL_NEIGHBORHOOD_SIZE)
print("CPU_THREADS:", CPU_THREADS)
if threadpool_info is not None:
    print("threadpool info:")
    for info in threadpool_info():
        print(" -", info.get("internal_api"), info.get("num_threads"), info.get("filepath", ""))

# %% [markdown]
# ## 2. 检查输入文件

# %%
# ==============================================================================
# 2. 检查输入文件
# ==============================================================================

def list_available_project_files(work_dir):
    """当关键输入缺失时，列出项目中可用的候选文件，避免编造路径。"""
    patterns = ["*.h5ad", "*.qs", "*.rds", "*.mtx", "*.txt", "*.csv"]
    found = []
    for pattern in patterns:
        found.extend(sorted(work_dir.rglob(pattern)))
    return found

missing = []
for p in [INPUT_H5AD, CILIAHUB_CSV]:
    if not p.exists():
        missing.append(p)

if missing:
    print("以下必需输入文件不存在：")
    for p in missing:
        print(" -", p)

    print("\n当前项目目录中发现的可用候选文件：")
    for p in list_available_project_files(WORK_DIR)[:200]:
        print(" -", p)

    raise FileNotFoundError(
        "请先确认上游 AnnData / CiliaHub CSV 路径。推荐输入是已有的 "
        "Malignant_RNA_assay_for_cNMF_countsX.h5ad。"
    )

print("输入文件检查通过。")
print("h5ad size:", INPUT_H5AD.stat().st_size / 1024**2, "MB")
print("CiliaHub CSV size:", CILIAHUB_CSV.stat().st_size / 1024, "KB")

# %% [markdown]
# ## 3. 读取已有 AnnData 对象

# %%
# ==============================================================================
# 3. 读取已有 AnnData 对象
# ==============================================================================

print("正在读取 AnnData...")
adata = ad.read_h5ad(INPUT_H5AD)

print("AnnData:", adata)
print("layers:", list(adata.layers.keys()))
print("obsm:", list(adata.obsm.keys()))
print("obs columns 前 30 个:", list(adata.obs.columns[:30]))
print("var columns:", list(adata.var.columns))

# cNMF 推荐使用 raw counts；当前项目导出的 h5ad 已把 raw counts 放在 X 和 layers['counts']
if "counts" in adata.layers:
    counts_all = adata.layers["counts"]
    print("使用 adata.layers['counts'] 作为 raw counts。")
else:
    counts_all = adata.X
    print("未检测到 layers['counts']，使用 adata.X 作为 counts。")

if not sp.issparse(counts_all):
    counts_all = sp.csr_matrix(counts_all)
else:
    counts_all = counts_all.tocsr()

print("counts shape cells x genes:", counts_all.shape)
print("counts dtype:", counts_all.dtype)
print("counts nnz:", counts_all.nnz)

# %% [markdown]
# ## 4. 读取 CiliaHub 基因并自动识别 gene symbol 列

# %%
# ==============================================================================
# 4. 读取 CiliaHub 基因并自动识别 gene symbol 列
# ==============================================================================

cilia_df = pd.read_csv(CILIAHUB_CSV, dtype=str)
print("CiliaHub CSV columns:", list(cilia_df.columns))
print(cilia_df.head())


def choose_gene_symbol_column(df):
    """自动选择最可能代表 gene symbol 的列。"""
    exact_priority = [
        "gene", "genes", "gene_symbol", "symbol", "hgnc_symbol",
        "external_gene_name", "gene.name", "gene_name"
    ]
    normalized = {c: re.sub(r"[^a-z0-9]+", "_", str(c).strip().lower()).strip("_") for c in df.columns}

    # 1) 优先匹配明确列名
    for target in exact_priority:
        for col, norm in normalized.items():
            if norm == target:
                return col, "列名精确匹配"

    # 2) 再按列名关键词和内容形态打分
    scores = []
    for col in df.columns:
        norm = normalized[col]
        values = df[col].dropna().astype(str).str.strip()
        values = values[values != ""]
        sample = values.head(200)
        if len(sample) == 0:
            scores.append((col, -999, "空列"))
            continue

        score = 0
        if "gene" in norm:
            score += 3
        if "symbol" in norm:
            score += 3
        if "ensembl" in norm or norm in {"ensg", "ensembl_id"}:
            score -= 4

        # HGNC symbol 通常为大写字母数字和少量符号，不像 ENSG ID，也不常包含逗号分隔多个 ID
        symbol_like = sample.str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$").mean()
        ensembl_like = sample.str.contains(r"ENSG\d+", regex=True).mean()
        comma_like = sample.str.contains(",").mean()
        score += 4 * symbol_like
        score -= 5 * ensembl_like
        score -= 2 * comma_like
        scores.append((col, score, f"symbol_like={symbol_like:.2f}, ensembl_like={ensembl_like:.2f}"))

    scores_df = pd.DataFrame(scores, columns=["column", "score", "note"]).sort_values("score", ascending=False)
    print("\n列名自动打分：")
    print(scores_df)

    if scores_df.empty or scores_df.iloc[0]["score"] < 3:
        raise ValueError(
            "无法可靠判断 gene symbol 列。请查看上面打印的列名，并手动设置 GENE_SYMBOL_COL。"
        )

    if len(scores_df) > 1 and scores_df.iloc[0]["score"] - scores_df.iloc[1]["score"] < 1:
        raise ValueError(
            "gene symbol 列判断不够明确。请查看上面打印的列名，并手动设置 GENE_SYMBOL_COL。"
        )

    return scores_df.iloc[0]["column"], "内容形态打分最高"


GENE_SYMBOL_COL, reason = choose_gene_symbol_column(cilia_df)
print(f"\n自动选择 gene symbol 列: {GENE_SYMBOL_COL} ({reason})")

ciliahub_genes = (
    cilia_df[GENE_SYMBOL_COL]
    .dropna()
    .astype(str)
    .str.strip()
)
ciliahub_genes = ciliahub_genes[ciliahub_genes != ""]
ciliahub_genes = pd.Index(pd.unique(ciliahub_genes), name="gene_symbol")

print("CiliaHub 去重后基因数:", len(ciliahub_genes))
print("前 20 个基因:", list(ciliahub_genes[:20]))

pd.DataFrame({"gene_symbol": ciliahub_genes}).to_csv(
    TABLES_DIR / "CiliaHub_gene_symbols_used.csv",
    index=False
)

# %% [markdown]
# ## 5. 检查 CiliaHub 基因与 AnnData var_names / gene_symbol 的交集

# %%
# ==============================================================================
# 5. 检查 CiliaHub 基因与 AnnData var_names / gene_symbol 的交集
# ==============================================================================

# 当前项目 h5ad 的 var_names 来自 gene_id；var['gene_symbol'] 保存原始 symbol。
# 因此同时检查 var_names 和 gene_symbol，优先用 gene_symbol 匹配，输出更容易解释。
var_names = pd.Index(adata.var_names.astype(str))

if "gene_symbol" in adata.var.columns:
    adata_gene_symbols = adata.var["gene_symbol"].astype(str)
    print("检测到 adata.var['gene_symbol']，优先用 gene_symbol 匹配。")
else:
    adata_gene_symbols = pd.Series(var_names.astype(str), index=adata.var_names)
    print("未检测到 adata.var['gene_symbol']，退回使用 adata.var_names 匹配。")

intersect_var_names = pd.Index(sorted(set(ciliahub_genes).intersection(set(var_names))))
intersect_gene_symbols = pd.Index(sorted(set(ciliahub_genes).intersection(set(adata_gene_symbols))))

print("CiliaHub 与 adata.var_names 交集基因数:", len(intersect_var_names))
print("CiliaHub 与 adata.var['gene_symbol'] 交集基因数:", len(intersect_gene_symbols))

if len(intersect_gene_symbols) >= len(intersect_var_names):
    keep_var_mask = adata_gene_symbols.isin(intersect_gene_symbols).values
    match_mode = "gene_symbol"
else:
    keep_var_mask = var_names.isin(intersect_var_names)
    match_mode = "var_names"

print("最终匹配模式:", match_mode)
print("最终保留基因数:", int(np.sum(keep_var_mask)))

if np.sum(keep_var_mask) < 20:
    missing_preview = sorted(set(ciliahub_genes) - set(adata_gene_symbols))[:30]
    raise ValueError(
        "CiliaHub 与表达矩阵交集少于 20 个基因，无法稳定运行 cNMF。\n"
        f"缺失基因示例: {missing_preview}\n"
        "请检查 CiliaHub CSV 使用的是 gene symbol 还是 Ensembl ID，以及 AnnData var 信息。"
    )

intersect_table = pd.DataFrame({
    "gene_symbol": adata_gene_symbols[keep_var_mask].values,
    "var_name": var_names[keep_var_mask].values,
})
intersect_table.to_csv(TABLES_DIR / "CiliaHub_intersect_genes_for_cNMF.csv", index=False)
intersect_table.head()

# %% [markdown]
# ## 6. 过滤低质量细胞和低表达 CiliaHub 基因

# %%
# ==============================================================================
# 6. 过滤低质量细胞和低表达 CiliaHub 基因
# ==============================================================================

# 先用全基因 raw counts 计算细胞 QC，避免 CiliaHub 限定基因集导致误判
cell_total_counts = np.asarray(counts_all.sum(axis=1)).ravel()
cell_detected_genes = np.asarray((counts_all > 0).sum(axis=1)).ravel()

cell_qc = pd.DataFrame({
    "total_counts_all_genes": cell_total_counts,
    "detected_genes_all_genes": cell_detected_genes,
}, index=adata.obs_names)
cell_qc.to_csv(TABLES_DIR / "cell_qc_all_genes.csv")

keep_cells = (
    (cell_total_counts >= MIN_CELL_TOTAL_COUNTS) &
    (cell_detected_genes >= MIN_CELL_DETECTED_GENES)
)

print("过滤前细胞数:", adata.n_obs)
print("保留细胞数:", int(np.sum(keep_cells)))
print("过滤掉细胞数:", int(np.sum(~keep_cells)))

# 再提取 CiliaHub 基因子集
adata_cilia = adata[keep_cells, keep_var_mask].copy()

# 用 raw counts 覆盖 X，确保 cNMF 输入是 count matrix
counts_cilia = counts_all[keep_cells, :][:, keep_var_mask].tocsr().astype(np.float32)
adata_cilia.X = counts_cilia
adata_cilia.layers["counts"] = counts_cilia.copy()

# 为了 cNMF 输出更好读，把 var_names 改成 gene symbol；若重复则 make_unique
if "gene_symbol" in adata_cilia.var.columns:
    adata_cilia.var["original_var_name"] = adata_cilia.var_names.astype(str)
    adata_cilia.var_names = pd.Index(adata_cilia.var["gene_symbol"].astype(str).values)
    adata_cilia.var_names_make_unique()

# 过滤低表达 CiliaHub 基因
counts_cilia = adata_cilia.layers["counts"].tocsr()
gene_total_counts = np.asarray(counts_cilia.sum(axis=0)).ravel()
gene_detected_cells = np.asarray((counts_cilia > 0).sum(axis=0)).ravel()

keep_genes = (
    (gene_total_counts >= MIN_GENE_TOTAL_COUNTS) &
    (gene_detected_cells >= MIN_GENE_EXPRESSED_CELLS)
)

gene_qc = pd.DataFrame({
    "gene": adata_cilia.var_names,
    "total_counts": gene_total_counts,
    "expressed_cells": gene_detected_cells,
    "keep": keep_genes,
})
gene_qc.to_csv(TABLES_DIR / "CiliaHub_gene_qc_for_cNMF.csv", index=False)

print("CiliaHub 交集基因数:", adata_cilia.n_vars)
print("低表达过滤后保留基因数:", int(np.sum(keep_genes)))

adata_cilia = adata_cilia[:, keep_genes].copy()
adata_cilia.layers["counts"] = adata_cilia.layers["counts"].tocsr().astype(np.float32)
adata_cilia.X = adata_cilia.layers["counts"].copy()

if adata_cilia.n_vars < max(20, K_VALUES.max() * 3):
    raise ValueError(
        f"过滤后仅剩 {adata_cilia.n_vars} 个基因，相对 K 最大值 {K_VALUES.max()} 太少。"
        "请降低 MIN_GENE_EXPRESSED_CELLS / MIN_GENE_TOTAL_COUNTS 或缩小 K_VALUES。"
    )

print("最终 cNMF 输入 AnnData:", adata_cilia)
print("最终基因示例:", list(adata_cilia.var_names[:20]))

# %% [markdown]
# ## 7. 保存 cNMF 输入矩阵

# %%
# ==============================================================================
# 7. 保存 cNMF 输入矩阵
# ==============================================================================

input_h5ad = OUT_DIR / "GSE132465_malignant_CiliaHub_counts_for_cNMF.h5ad"
counts_txt = OUT_DIR / "GSE132465_malignant_CiliaHub_counts_for_cNMF.txt"

print("正在保存 h5ad:", input_h5ad)
adata_cilia.write_h5ad(input_h5ad, compression="gzip")

# cNMF prepare 最稳定的输入格式是 cells x genes 的 tab-delimited count matrix
# CiliaHub 子集规模较小，转成 dense DataFrame 便于 cNMF 读取和人工检查
X = adata_cilia.layers["counts"]
if sp.issparse(X):
    X = X.toarray()

counts_df = pd.DataFrame(
    X,
    index=adata_cilia.obs_names.astype(str),
    columns=adata_cilia.var_names.astype(str),
)

# 避免 -0.0；counts 应为非负
counts_df[counts_df < 0] = 0

print("正在保存 counts txt:", counts_txt)
counts_df.to_csv(counts_txt, sep="\t")

print("counts_df shape:", counts_df.shape)
print("输出文件：")
print(" -", input_h5ad)
print(" -", counts_txt)

# %% [markdown]
# ## 8. 运行 cNMF：prepare / factorize / combine / K selection

# %%
# ==============================================================================
# 8. 运行 cNMF
# ==============================================================================

# 这一单元耗时最长。N_ITER=100 且 K=2..12 时通常需要较久。
# 调试时可以先把 N_ITER 改成 20，确认流程后再改回 100 或更高。

cnmf_obj = cNMF(output_dir=str(OUT_DIR), name=CNMF_NAME)

print("开始 cNMF prepare...")
try:
    # 限定到 CiliaHub 后，尽量让 cNMF 使用过滤后的全部基因，而不是再默认筛 2000 个高变基因。
    cnmf_obj.prepare(
        counts_fn=str(counts_txt),
        components=K_VALUES,
        n_iter=N_ITER,
        seed=RANDOM_SEED,
        num_highvar_genes=adata_cilia.n_vars,
    )
except TypeError:
    # 兼容旧版本 cnmf：如果没有 num_highvar_genes 参数，则退回基础参数。
    cnmf_obj.prepare(
        counts_fn=str(counts_txt),
        components=K_VALUES,
        n_iter=N_ITER,
        seed=RANDOM_SEED,
    )

print("开始 cNMF factorize...")
cnmf_obj.factorize(worker_i=WORKER_INDEX, total_workers=TOTAL_WORKERS, skip_completed_runs=True)

print("开始 cNMF combine...")
try:
    cnmf_obj.combine(skip_missing_files=False)
except TypeError:
    cnmf_obj.combine()

print("生成 K selection plot...")
cnmf_obj.k_selection_plot()

print("cNMF K selection 完成。请查看：")
print(CNMF_RUN_DIR)
print("常见图文件：", CNMF_RUN_DIR / f"{CNMF_NAME}.k_selection.png")

# %%
print(adata.shape)

# %% [markdown]
# ## 9. 查看 K selection 结果并设定最终 K

# %%
# ==============================================================================
# 9. 查看 K selection 结果并设定最终 K
# ==============================================================================

# cNMF 不同版本的 K selection 统计文件名可能略有差异，这里自动搜索并打印。
k_selection_files = sorted(CNMF_RUN_DIR.glob("*k_selection*"))
print("K selection 相关文件：")
for p in k_selection_files:
    print(" -", p.name)

# 如果 SELECTED_K 仍为 None，默认取 K_VALUES 中间值，便于 notebook 可继续运行。
# 实际分析建议先查看 k_selection.png，再手动把 SELECTED_K 改成拐点或稳定性合适的 K。
if SELECTED_K is None:
    selected_k = int(K_VALUES[len(K_VALUES) // 2])
    print(f"SELECTED_K=None，临时使用中间 K={selected_k}。正式结果请结合 K selection plot 手动确认。")
else:
    selected_k = int(SELECTED_K)
    print("使用手动指定 SELECTED_K:", selected_k)

if selected_k not in set(map(int, K_VALUES)):
    raise ValueError(f"selected_k={selected_k} 不在 K_VALUES={list(K_VALUES)} 中。")

# %% [markdown]
# ## 10. 运行 consensus 并保存核心结果

# %%
# ==============================================================================
# 10. 运行 consensus 并保存核心结果
# ==============================================================================

cnmf_obj = cNMF(output_dir=str(OUT_DIR), name=CNMF_NAME)

# 如果之前用过不同 local_neighborhood_size 跑失败，旧 local density cache 会影响重跑；这里删除当前 K 的缓存。
for density_cache in CNMF_RUN_DIR.glob(f"cnmf_tmp/*local_density_cache.k_{selected_k}*.npz"):
    print("删除旧 density cache:", density_cache)
    density_cache.unlink()

print(
    f"开始 consensus: K={selected_k}, density_threshold={DENSITY_THRESHOLD}, "
    f"local_neighborhood_size={LOCAL_NEIGHBORHOOD_SIZE}"
)
try:
    consensus_result = cnmf_obj.consensus(
        k=selected_k,
        density_threshold=DENSITY_THRESHOLD,
        local_neighborhood_size=LOCAL_NEIGHBORHOOD_SIZE,
        show_clustering=True,
    )
except TypeError:
    # 兼容部分版本参数名
    consensus_result = cnmf_obj.consensus(
        k=selected_k,
        local_density_threshold=DENSITY_THRESHOLD,
        local_neighborhood_size=LOCAL_NEIGHBORHOOD_SIZE,
        show_clustering=True,
    )

usage = spectra_scores = spectra_tpm = top_genes = None
if isinstance(consensus_result, tuple):
    if len(consensus_result) >= 1:
        usage = consensus_result[0]
    if len(consensus_result) >= 2:
        spectra_scores = consensus_result[1]
    if len(consensus_result) >= 3:
        spectra_tpm = consensus_result[2]
    if len(consensus_result) >= 4:
        top_genes = consensus_result[3]

# cNMF 输出文件名中常把 0.01 写成 dt_0_01；这里用通配符自动收集，避免版本差异。
consensus_files = sorted(CNMF_RUN_DIR.glob(f"*k_{selected_k}*")) + sorted(CNMF_RUN_DIR.glob(f"*k{selected_k}*"))
consensus_files = sorted(set(consensus_files))

print("consensus 相关输出文件：")
for p in consensus_files:
    print(" -", p.name)

# 如果当前版本 consensus 返回了对象，额外保存一份到 tables 目录，方便后续读取
if usage is not None:
    usage_df = pd.DataFrame(usage)
    usage_df.to_csv(TABLES_DIR / f"usage_matrix.k_{selected_k}.csv")
    print("已保存 usage:", TABLES_DIR / f"usage_matrix.k_{selected_k}.csv")

if spectra_scores is not None:
    spectra_scores_df = pd.DataFrame(spectra_scores)
    spectra_scores_df.to_csv(TABLES_DIR / f"consensus_programs_gene_spectra_score.k_{selected_k}.csv")
    print("已保存 spectra_scores:", TABLES_DIR / f"consensus_programs_gene_spectra_score.k_{selected_k}.csv")

if spectra_tpm is not None:
    spectra_tpm_df = pd.DataFrame(spectra_tpm)
    spectra_tpm_df.to_csv(TABLES_DIR / f"consensus_programs_gene_spectra_tpm.k_{selected_k}.csv")
    print("已保存 spectra_tpm:", TABLES_DIR / f"consensus_programs_gene_spectra_tpm.k_{selected_k}.csv")

if top_genes is not None:
    pd.DataFrame(top_genes).to_csv(TABLES_DIR / f"top_genes_raw_from_cNMF.k_{selected_k}.csv")
    print("已保存 top_genes:", TABLES_DIR / f"top_genes_raw_from_cNMF.k_{selected_k}.csv")

# %% [markdown]
# ## 11. 读取 cNMF 输出并生成 marker gene 排名

# %%
# ==============================================================================
# 11. 读取 cNMF 输出并生成 marker gene 排名
# ==============================================================================


def find_one(patterns):
    """按多个通配符寻找一个结果文件。"""
    hits = []
    for pat in patterns:
        hits.extend(CNMF_RUN_DIR.glob(pat))
    hits = sorted(set(hits))
    if len(hits) == 0:
        return None
    return hits[-1]

usage_file = find_one([f"*.usages.k_{selected_k}*.txt", f"*.usages.k{selected_k}*.txt"])
spectra_score_file = find_one([f"*.gene_spectra_score.k_{selected_k}*.txt", f"*.gene_spectra_score.k{selected_k}*.txt"])
spectra_tpm_file = find_one([f"*.gene_spectra_tpm.k_{selected_k}*.txt", f"*.gene_spectra_tpm.k{selected_k}*.txt"])
top_genes_file = find_one([f"*.top_genes.k_{selected_k}*.txt", f"*.top_genes.k{selected_k}*.txt"])

print("usage_file:", usage_file)
print("spectra_score_file:", spectra_score_file)
print("spectra_tpm_file:", spectra_tpm_file)
print("top_genes_file:", top_genes_file)

if usage_file is None or spectra_score_file is None:
    raise FileNotFoundError(
        "没有找到 usage 或 gene_spectra_score 文件。请确认 consensus 单元是否成功运行，"
        "并检查 CNMF_RUN_DIR 下的实际文件名。"
    )

usage_df = pd.read_csv(usage_file, sep="\t", index_col=0)
spectra_score_df = pd.read_csv(spectra_score_file, sep="\t", index_col=0)

if spectra_tpm_file is not None:
    spectra_tpm_df = pd.read_csv(spectra_tpm_file, sep="\t", index_col=0)
else:
    spectra_tpm_df = None

print("usage_df shape:", usage_df.shape)
print("spectra_score_df shape:", spectra_score_df.shape)
print("usage columns:", list(usage_df.columns))
print("spectra_score index 前几个:", list(spectra_score_df.index[:5]))
print("spectra_score columns 前几个:", list(spectra_score_df.columns[:5]))

# cNMF spectra 文件在不同版本中可能是 program x gene 或 gene x program。
# 这里根据与输入基因名的重合度自动判断方向，并统一成 program x gene。
input_genes = set(adata_cilia.var_names.astype(str))
index_gene_overlap = len(set(map(str, spectra_score_df.index)).intersection(input_genes))
column_gene_overlap = len(set(map(str, spectra_score_df.columns)).intersection(input_genes))

if column_gene_overlap >= index_gene_overlap:
    spectra_program_by_gene = spectra_score_df.copy()
else:
    spectra_program_by_gene = spectra_score_df.T.copy()

spectra_program_by_gene.index = [f"Program_{i+1}" for i in range(spectra_program_by_gene.shape[0])]

marker_rows = []
for program in spectra_program_by_gene.index:
    ranked = spectra_program_by_gene.loc[program].sort_values(ascending=False)
    for rank, (gene, score) in enumerate(ranked.items(), start=1):
        marker_rows.append({
            "program": program,
            "rank": rank,
            "gene": gene,
            "spectra_score": score,
        })

marker_rank_df = pd.DataFrame(marker_rows)
marker_rank_path = TABLES_DIR / f"marker_gene_rankings_by_program.k_{selected_k}.csv"
marker_rank_df.to_csv(marker_rank_path, index=False)

# 同步保存标准化命名后的核心矩阵，方便后续分析
usage_df.to_csv(TABLES_DIR / f"usage_matrix_from_cNMF.k_{selected_k}.csv")
spectra_program_by_gene.to_csv(TABLES_DIR / f"consensus_programs_by_gene.k_{selected_k}.csv")

print("marker gene 排名已保存:", marker_rank_path)
marker_rank_df.head(20)

# %% [markdown]
# ## 12. 基础可视化：usage heatmap 和 program gene heatmap

# %%
# ==============================================================================
# 12. 基础可视化
# ==============================================================================

if sns is None:
    print("未安装 seaborn，跳过 heatmap。")
else:
    # --------------------------------------------------------------------------
    # 12.1 单细胞 usage heatmap
    # --------------------------------------------------------------------------
    usage_plot_df = usage_df.copy()
    usage_plot_df.columns = [f"Program_{i+1}" for i in range(usage_plot_df.shape[1])]

    # 细胞很多时，完整 clustermap 会较慢；这里按 usage 主程序排序，并最多抽样 3000 个细胞展示
    max_cells_plot = 3000
    usage_order = usage_plot_df.idxmax(axis=1).sort_values().index
    usage_plot_df = usage_plot_df.loc[usage_order]
    if usage_plot_df.shape[0] > max_cells_plot:
        rng = np.random.default_rng(RANDOM_SEED)
        sampled = rng.choice(usage_plot_df.index.values, size=max_cells_plot, replace=False)
        usage_plot_df = usage_plot_df.loc[sampled].sort_index()

    plt.figure(figsize=(max(6, selected_k * 0.7), 8))
    sns.heatmap(
        usage_plot_df,
        cmap="viridis",
        xticklabels=True,
        yticklabels=False,
    )
    plt.title(f"cNMF usage heatmap, K={selected_k}")
    plt.xlabel("Program")
    plt.ylabel("Cells")
    plt.tight_layout()
    usage_heatmap_path = PLOTS_DIR / f"usage_heatmap.k_{selected_k}.png"
    plt.savefig(usage_heatmap_path, dpi=300)
    plt.show()
    print("已保存:", usage_heatmap_path)

    # --------------------------------------------------------------------------
    # 12.2 每个 program 的 top genes spectra heatmap
    # --------------------------------------------------------------------------
    top_n = 15
    top_genes_union = []
    for program in spectra_program_by_gene.index:
        top_genes_union.extend(
            spectra_program_by_gene.loc[program].sort_values(ascending=False).head(top_n).index.tolist()
        )
    top_genes_union = list(dict.fromkeys(top_genes_union))

    program_gene_mat = spectra_program_by_gene.loc[:, top_genes_union]

    plt.figure(figsize=(max(10, len(top_genes_union) * 0.22), max(3, selected_k * 0.45)))
    sns.heatmap(
        program_gene_mat,
        cmap="mako",
        xticklabels=True,
        yticklabels=True,
    )
    plt.title(f"Top CiliaHub genes per cNMF program, K={selected_k}")
    plt.xlabel("Genes")
    plt.ylabel("Programs")
    plt.tight_layout()
    program_heatmap_path = PLOTS_DIR / f"program_top_gene_heatmap.k_{selected_k}.png"
    plt.savefig(program_heatmap_path, dpi=300)
    plt.show()
    print("已保存:", program_heatmap_path)

# %% [markdown]
# ## 13. 可选：按已有 metadata 汇总 program usage

# %%
# ==============================================================================
# 13. 可选：按已有 metadata 汇总 program usage
# ==============================================================================

# 自动寻找常见分组列；当前对象通常有 seurat_clusters / Patient / orig.ident 等列
candidate_group_cols = [
    "seurat_clusters", "Patient", "orig.ident", "Class", "celltype", "cell_type", "sample", "Sample"
]
group_cols = [c for c in candidate_group_cols if c in adata_cilia.obs.columns]
print("可用于汇总 usage 的 metadata 列:", group_cols)

if len(group_cols) == 0:
    print("没有找到合适的分组列，跳过分组 usage heatmap。")
elif sns is None:
    print("未安装 seaborn，跳过分组 usage heatmap。")
else:
    group_col = group_cols[0]
    print("默认使用分组列:", group_col)

    usage_plot_df = usage_df.copy()
    usage_plot_df.columns = [f"Program_{i+1}" for i in range(usage_plot_df.shape[1])]

    obs_group = adata_cilia.obs.loc[usage_plot_df.index, group_col].astype(str)
    usage_by_group = usage_plot_df.groupby(obs_group).mean()

    usage_by_group.to_csv(TABLES_DIR / f"usage_mean_by_{group_col}.k_{selected_k}.csv")

    plt.figure(figsize=(max(6, selected_k * 0.7), max(4, usage_by_group.shape[0] * 0.35)))
    sns.heatmap(
        usage_by_group,
        cmap="viridis",
        xticklabels=True,
        yticklabels=True,
    )
    plt.title(f"Mean cNMF usage by {group_col}, K={selected_k}")
    plt.xlabel("Program")
    plt.ylabel(group_col)
    plt.tight_layout()
    group_usage_path = PLOTS_DIR / f"usage_mean_by_{group_col}.k_{selected_k}.png"
    plt.savefig(group_usage_path, dpi=300)
    plt.show()
    print("已保存:", group_usage_path)

# %% [markdown]
# ## 14. 输出文件汇总

# %%
# ==============================================================================
# 14. 输出文件汇总
# ==============================================================================

print("主要输出目录:")
print("OUT_DIR:", OUT_DIR)
print("CNMF_RUN_DIR:", CNMF_RUN_DIR)
print("TABLES_DIR:", TABLES_DIR)
print("PLOTS_DIR:", PLOTS_DIR)

print("\nOUT_DIR 文件：")
for p in sorted(OUT_DIR.glob("*")):
    print(" -", p)

print("\nTABLES_DIR 文件：")
for p in sorted(TABLES_DIR.glob("*")):
    print(" -", p.name)

print("\nPLOTS_DIR 文件：")
for p in sorted(PLOTS_DIR.glob("*")):
    print(" -", p.name)

# %% [markdown]
# ## 15. 运行顺序和依赖说明
# 
# 推荐运行顺序：
# 
# 1. 先运行第 1-7 节，确认输入对象、CiliaHub gene symbol 列、基因交集、过滤后细胞数和基因数。
# 2. 运行第 8 节执行 cNMF。正式结果默认 `N_ITER = 100`、`K_VALUES = 2..12`；如果只是检查流程，可在启动前设置 `CNMF_FAST_DEBUG=1`，会写入单独的 `cNMF_CiliaHub_GSE132465_debug/` 目录。
# 3. 查看 `cNMF_CiliaHub_GSE132465/GSE132465_malignant_CiliaHub_cNMF/GSE132465_malignant_CiliaHub_cNMF.k_selection.png`，然后在第 1 节把 `SELECTED_K` 改成最终 K。
# 4. 运行第 10-14 节生成 consensus programs、usage matrix、marker gene ranking 和基础可视化。
# 
# 需要的 Python 包：
# 
# ```bash
# pip install cnmf anndata pandas scipy scikit-learn scanpy seaborn matplotlib
# ```
# 
# 如果用 conda，更推荐单独建环境，例如：
# 
# ```bash
# conda create -n cnmf_env -c conda-forge -c bioconda python=3.10 cnmf scanpy anndata pandas scipy scikit-learn seaborn matplotlib
# conda activate cnmf_env
# ```
# 
# 本项目当前已在 `CRCpy` conda 环境中完成快速调试运行。命令示例：
# 
# ```bash
# CNMF_FAST_DEBUG=1 NUMBA_CACHE_DIR=/tmp/numba-cache MPLCONFIGDIR=/tmp/matplotlib-cnmf /home/zerui/miniconda3/envs/CRCpy/bin/python
# ```
# 
# 主要输出文件：
# 
# - `GSE132465_malignant_CiliaHub_counts_for_cNMF.h5ad`：过滤后的 cNMF 输入 AnnData。
# - `GSE132465_malignant_CiliaHub_counts_for_cNMF.txt`：cells x genes 的 raw count matrix，供 cNMF 使用。
# - `tables/CiliaHub_intersect_genes_for_cNMF.csv`：CiliaHub 与表达矩阵交集基因。
# - `tables/CiliaHub_gene_qc_for_cNMF.csv`：CiliaHub 基因过滤统计。
# - `tables/usage_matrix_from_cNMF.k_*.csv`：细胞 x program 的 usage matrix。
# - `tables/consensus_programs_by_gene.k_*.csv`：program x gene 的 consensus program 矩阵。
# - `tables/marker_gene_rankings_by_program.k_*.csv`：每个 program 的 marker gene 排名。
# - `plots/usage_heatmap.k_*.png`：usage heatmap。
# - `plots/program_top_gene_heatmap.k_*.png`：program top gene heatmap。
# - `GSE132465_malignant_CiliaHub_cNMF/*.k_selection.png`：cNMF 自动生成的 K 选择图。


