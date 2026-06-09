#!/usr/bin/env python3
"""
09_scvelo_velocity_combined_introns.py
Run scVelo RNA velocity for combined Aorta Macrophage + PBMC Monocyte cells.

Inputs are produced by 08_combined_macrophage_monocyte_monocle2_velocity_introns.R:
  files/Combined/AortaMacrophage_PBMCMonocyte/*_velocity_metadata.csv
  files/Combined/AortaMacrophage_PBMCMonocyte/*_loom_file_map.csv

Use:
  /home/zerui/miniconda3/envs/velocity/bin/python \
    /home/zerui/code/coldmouse_introns/scripts/09_scvelo_velocity_combined_introns.py
"""

import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba_cache_velocity")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplconfig_velocity")
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)

import gc
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as plt

ROOT = Path(os.environ.get("COLDMOUSE_INTRONS_OUTPUT", "/mnt/disk1/qiuzerui/expriments/coldmouse_introns"))
TAG = "AortaMacrophage_PBMCMonocyte"
FILES_DIR = ROOT / "files" / "Combined" / TAG
PLOTS_DIR = ROOT / "plots" / "Combined" / TAG
RAWDATA_DIR = ROOT / "rawdata"
FILES_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

METADATA_FILE = Path(os.environ.get(
    "SCVELO_METADATA",
    FILES_DIR / f"introns_combined_{TAG}_velocity_metadata.csv",
))
LOOM_MAP_FILE = Path(os.environ.get(
    "SCVELO_LOOM_MAP",
    FILES_DIR / f"introns_combined_{TAG}_loom_file_map.csv",
))
MIN_SHARED_COUNTS = int(os.environ.get("SCVELO_MIN_SHARED_COUNTS", "20"))
N_TOP_GENES = int(os.environ.get("SCVELO_N_TOP_GENES", "2000"))
N_PCS = int(os.environ.get("SCVELO_N_PCS", "30"))
N_NEIGHBORS = int(os.environ.get("SCVELO_N_NEIGHBORS", "30"))
VELOCITY_MODE = os.environ.get("SCVELO_MODE", "stochastic")


def read_subset_loom(loom_file: Path, wanted_obs: set):
    print(f"Reading loom: {loom_file}")
    ad = sc.read_loom(str(loom_file), sparse=True)
    if "obs_names" in ad.obs.columns:
        ad.obs_names = ad.obs["obs_names"].astype(str).values
    overlap = ad.obs_names.intersection(list(wanted_obs))
    print(f"  cells in metadata: {len(wanted_obs)}; matched in loom: {len(overlap)}")
    if len(overlap) == 0:
        return None
    ad = ad[overlap, :].copy()
    ad.obs["loom_file"] = loom_file.name
    return ad


def save_velocity_plot(adata, color, filename, stream=True):
    plt.figure()
    if stream:
        scv.pl.velocity_embedding_stream(
            adata,
            basis="umap",
            color=color,
            legend_loc="right margin",
            title=f"scVelo {color}",
            show=False,
        )
    else:
        scv.pl.velocity_embedding(
            adata,
            basis="umap",
            color=color,
            arrow_length=3,
            arrow_size=2,
            dpi=120,
            show=False,
        )
    out_png = PLOTS_DIR / f"{filename}.png"
    out_pdf = PLOTS_DIR / f"{filename}.pdf"
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    try:
        plt.savefig(out_pdf, bbox_inches="tight")
    except ValueError as exc:
        print(f"WARNING: vector PDF save failed for {out_pdf}: {exc}")
        plt.close("all")
        img = plt.imread(out_png)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(img)
        ax.axis("off")
        fig.savefig(out_pdf, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved raster PDF fallback: {out_pdf}")
    else:
        plt.close("all")
    print(f"Saved: {out_png}")


def main():
    print("scVelo version:", scv.__version__)
    print("metadata:", METADATA_FILE)
    print("loom map:", LOOM_MAP_FILE)
    if not METADATA_FILE.exists():
        raise FileNotFoundError(f"metadata not found: {METADATA_FILE}. Run script 08 first.")
    if not LOOM_MAP_FILE.exists():
        raise FileNotFoundError(f"loom map not found: {LOOM_MAP_FILE}. Run script 08 first.")

    meta = pd.read_csv(METADATA_FILE)
    required = {"loom_obs_name", "seurat_cell", "umap_1", "umap_2", "sample_id"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"metadata missing columns: {sorted(missing)}")
    meta = meta.drop_duplicates("loom_obs_name").set_index("loom_obs_name", drop=False)
    loom_map = pd.read_csv(LOOM_MAP_FILE)

    adatas = []
    for _, row in loom_map.iterrows():
        sample_id = str(row["sample_id"])
        loom_file = Path(row.get("loom_file", RAWDATA_DIR / f"{sample_id}.loom"))
        if not loom_file.exists():
            loom_file = RAWDATA_DIR / f"{sample_id}.loom"
        wanted = set(meta.loc[meta["sample_id"].astype(str) == sample_id, "loom_obs_name"])
        if not wanted:
            print(f"Skip {sample_id}: no cells in metadata")
            continue
        ad = read_subset_loom(loom_file, wanted)
        if ad is not None:
            adatas.append(ad)
        del ad
        gc.collect()

    if not adatas:
        raise RuntimeError("No loom cells matched metadata. Check loom_obs_name format.")

    print("Concatenating loom subsets...")
    adata = adatas[0].concatenate(
        adatas[1:],
        join="outer",
        batch_key="loom_batch",
        index_unique=None,
    ) if len(adatas) > 1 else adatas[0]
    del adatas
    gc.collect()

    matched_meta = meta.loc[adata.obs_names].copy()
    for col in matched_meta.columns:
        if col in {"umap_1", "umap_2"}:
            continue
        adata.obs[col] = matched_meta[col].astype(str).values
    adata.obsm["X_umap"] = matched_meta[["umap_1", "umap_2"]].to_numpy(dtype=float)
    adata.obs_names = matched_meta["seurat_cell"].astype(str).values
    adata.obs_names_make_unique()

    print(adata)
    print("layers:", list(adata.layers.keys()))
    if "spliced" not in adata.layers or "unspliced" not in adata.layers:
        raise ValueError("spliced/unspliced layers not found in loom data.")

    scv.settings.verbosity = 3
    scv.settings.presenter_view = False
    scv.settings.set_figure_params("scvelo", dpi=120, dpi_save=300, transparent=False)

    print("filter_and_normalize...")
    scv.pp.filter_and_normalize(
        adata,
        min_shared_counts=MIN_SHARED_COUNTS,
    )
    print("moments...")
    scv.pp.moments(adata, n_pcs=N_PCS, n_neighbors=N_NEIGHBORS)
    print(f"velocity mode={VELOCITY_MODE}...")
    scv.tl.velocity(adata, mode=VELOCITY_MODE)
    scv.tl.velocity_graph(adata, n_jobs=1, backend="threading", show_progress_bar=False)
    scv.tl.velocity_confidence(adata)

    velocity_meta = adata.obs.copy()
    keep_cols = [c for c in [
        "seurat_cell", "loom_obs_name", "sample_id", "source_subset", "source_tissue",
        "source_celltype", "Group", "Condition", "combined_clusters",
        "velocity_length", "velocity_confidence",
    ] if c in velocity_meta.columns]
    velocity_meta = velocity_meta[keep_cols]
    velocity_meta.to_csv(FILES_DIR / f"introns_combined_{TAG}_scvelo_velocity_metadata.csv")

    h5ad_file = FILES_DIR / f"introns_combined_{TAG}_scvelo.h5ad"
    adata.write_h5ad(h5ad_file)
    print(f"Saved h5ad: {h5ad_file}")

    for key in ["X_umap", "velocity_umap"]:
        if key in adata.obsm:
            arr = np.asarray(adata.obsm[key], dtype=float)
            arr[~np.isfinite(arr)] = np.nan
            adata.obsm[key] = arr

    for color in ["source_subset", "Group", "combined_clusters"]:
        if color in adata.obs.columns:
            save_velocity_plot(adata, color, f"introns_combined_{TAG}_scvelo_stream_{color}", stream=True)
            save_velocity_plot(adata, color, f"introns_combined_{TAG}_scvelo_arrows_{color}", stream=False)

    print("scVelo analysis complete.")


if __name__ == "__main__":
    main()
