# %%
# ==============================================================================
# 02_make_h5ad_from_mtx.py
# 从 R 导出的 Matrix Market / TSV 文件构建 AnnData，并保存为 h5ad
# ==============================================================================

from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad


# ------------------------------------------------------------------------------
# 1. 路径设置
# ------------------------------------------------------------------------------

work_dir = Path("/mnt/disk1/qiuzerui/downloads/CRC/GSE132465")

in_dir = work_dir / "Malignant_RNA_assay_export_for_python"

out_h5ad = work_dir / "Malignant_RNA_assay_for_cNMF_countsX.h5ad"

print("输入目录:")
print(in_dir)

print("\n输出 h5ad:")
print(out_h5ad)


# ------------------------------------------------------------------------------
# 2. 读取 counts matrix
# ------------------------------------------------------------------------------

mtx_path = in_dir / "counts_genes_by_cells.mtx"

if not mtx_path.exists():
    raise FileNotFoundError(f"找不到 counts matrix: {mtx_path}")

print("\n正在读取 counts matrix...")

counts = scipy.io.mmread(mtx_path)

# R 导出的矩阵是 genes × cells
# AnnData 需要 cells × genes
if not sp.issparse(counts):
    counts = sp.coo_matrix(counts)

counts = counts.T.tocsr()

# 用 float32 节省空间；数值仍然是 raw counts
counts = counts.astype(np.float32)

print("counts shape cells × genes:", counts.shape)
print("counts nnz:", counts.nnz)


# ------------------------------------------------------------------------------
# 3. 读取 genes / barcodes / metadata
# ------------------------------------------------------------------------------

genes_path = in_dir / "genes.tsv"
barcodes_path = in_dir / "barcodes.tsv"
metadata_path = in_dir / "metadata.tsv"

for p in [genes_path, barcodes_path, metadata_path]:
    if not p.exists():
        raise FileNotFoundError(f"找不到文件: {p}")

print("\n正在读取 genes / barcodes / metadata...")

genes = pd.read_csv(genes_path, sep="\t", dtype=str)
barcodes = pd.read_csv(barcodes_path, sep="\t", dtype=str)
metadata = pd.read_csv(metadata_path, sep="\t", low_memory=False)

required_gene_cols = {"gene_id", "gene_symbol"}
if not required_gene_cols.issubset(set(genes.columns)):
    raise ValueError(f"genes.tsv 必须包含列: {required_gene_cols}")

if "cell" not in barcodes.columns:
    raise ValueError("barcodes.tsv 必须包含 cell 列")

if "cell" not in metadata.columns:
    raise ValueError("metadata.tsv 必须包含 cell 列")

cells = barcodes["cell"].astype(str).values

if counts.shape[0] != len(cells):
    raise ValueError(
        f"细胞数不匹配: counts 有 {counts.shape[0]} 个细胞, "
        f"barcodes.tsv 有 {len(cells)} 个细胞"
    )

if counts.shape[1] != genes.shape[0]:
    raise ValueError(
        f"基因数不匹配: counts 有 {counts.shape[1]} 个基因, "
        f"genes.tsv 有 {genes.shape[0]} 个基因"
    )


# ------------------------------------------------------------------------------
# 4. 整理 obs
# ------------------------------------------------------------------------------

print("\n正在整理 obs metadata...")

metadata["cell"] = metadata["cell"].astype(str)
metadata = metadata.set_index("cell")

missing_cells = set(cells) - set(metadata.index)
if len(missing_cells) > 0:
    example_missing = list(missing_cells)[:5]
    raise ValueError(f"metadata 中缺少部分细胞，例如: {example_missing}")

obs = metadata.loc[cells].copy()

# h5ad 对 object 类型比较敏感，统一清理字符列
for col in obs.columns:
    if obs[col].dtype == "object" or pd.api.types.is_string_dtype(obs[col]):
        obs[col] = obs[col].fillna("NA").astype(str)
    elif pd.api.types.is_bool_dtype(obs[col]):
        obs[col] = obs[col].fillna(False).astype(bool)

obs.index = pd.Index(cells, name=None)


# ------------------------------------------------------------------------------
# 5. 整理 var
# ------------------------------------------------------------------------------

print("正在整理 var gene metadata...")

genes["gene_id"] = genes["gene_id"].astype(str)
genes["gene_symbol"] = genes["gene_symbol"].astype(str)

var = genes.copy()
var.index = pd.Index(var["gene_id"].values, name=None)

# 保证 var index 唯一
if not var.index.is_unique:
    print("警告：gene_id 仍有重复，将自动 make_unique。")
    var.index = pd.Index(ad.utils.make_index_unique(var.index))


# ------------------------------------------------------------------------------
# 6. 构建 AnnData
# ------------------------------------------------------------------------------

print("\n正在构建 AnnData...")

adata = ad.AnnData(
    X=counts,
    obs=obs,
    var=var
)

adata.obs_names = pd.Index(cells, name=None)
adata.var_names = pd.Index(var.index.astype(str), name=None)

adata.var_names_make_unique()

# raw counts 同时放入 layers["counts"]
adata.layers["counts"] = adata.X.copy()

adata.uns["source"] = "Malignant_RNA_assay.qs"
adata.uns["X_matrix"] = "raw counts"
adata.uns["note"] = "Converted by R mtx export + Python AnnData construction"


# ------------------------------------------------------------------------------
# 7. 可选：读取 logdata
# ------------------------------------------------------------------------------

logdata_path = in_dir / "logdata_genes_by_cells.mtx"

if logdata_path.exists():
    print("\n检测到 logdata matrix，正在读取...")
    
    logdata = scipy.io.mmread(logdata_path)
    
    if not sp.issparse(logdata):
        logdata = sp.coo_matrix(logdata)
    
    logdata = logdata.T.tocsr().astype(np.float32)
    
    if logdata.shape == adata.shape:
        adata.layers["logdata"] = logdata
        print("已加入 adata.layers['logdata']")
    else:
        print("logdata shape 与 adata 不匹配，跳过。")


# ------------------------------------------------------------------------------
# 8. 读取 reductions
# ------------------------------------------------------------------------------

print("\n正在读取 reductions...")

for red_file in sorted(in_dir.glob("reduction_*.tsv")):
    red_name = red_file.stem.replace("reduction_", "")
    
    emb = pd.read_csv(red_file, sep="\t")
    
    if "cell" not in emb.columns:
        print(f"跳过 {red_file.name}: 没有 cell 列")
        continue
    
    emb["cell"] = emb["cell"].astype(str)
    emb = emb.set_index("cell")
    
    missing_cells = set(adata.obs_names) - set(emb.index)
    if len(missing_cells) > 0:
        print(f"跳过 {red_file.name}: 有细胞缺失")
        continue
    
    emb = emb.loc[adata.obs_names]
    
    # 所有坐标列转为数值
    emb_numeric = emb.apply(pd.to_numeric, errors="coerce")
    
    if emb_numeric.isna().any().any():
        print(f"警告：{red_file.name} 中存在 NA 坐标，将跳过该 reduction")
        continue
    
    key = f"X_{red_name}"
    adata.obsm[key] = emb_numeric.values.astype(np.float32)
    
    print(f"已加入 adata.obsm['{key}']: {emb_numeric.shape}")


# ------------------------------------------------------------------------------
# 9. 最终检查
# ------------------------------------------------------------------------------

print("\n================ AnnData 检查 ================")
print(adata)

print("\nobs 前几列:")
print(adata.obs.head())

print("\nvar 前几列:")
print(adata.var.head())

print("\nobsm keys:")
print(list(adata.obsm.keys()))

print("\nlayers:")
print(list(adata.layers.keys()))

print("\nX 是否稀疏矩阵:", sp.issparse(adata.X))
print("X shape:", adata.X.shape)
print("X dtype:", adata.X.dtype)


# ------------------------------------------------------------------------------
# 10. 保存 h5ad
# ------------------------------------------------------------------------------

print("\n正在写出 h5ad...")

adata.write_h5ad(out_h5ad, compression="gzip")

print("\n================ 转换完成 ================")
print("输出文件:")
print(out_h5ad)
