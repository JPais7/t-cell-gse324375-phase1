#!/usr/bin/env python3
"""Audit all filtered GEX matrices without applying cell filters."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scanpy as sc


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "extracted"
OUT = ROOT / "results" / "qc"
OUT.mkdir(parents=True, exist_ok=True)


def library_key(path: Path) -> str:
    name = path.name.split("_", 1)[1]
    return re.sub(r"_GEX(?:_FL)?_filtered_feature_bc_matrix\.h5$", "", name)


def matching_csv(gex: Path, modality: str) -> Path:
    key = library_key(gex)
    matches = list(DATA.glob(f"*_{key}_{modality}*.csv.gz"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {modality} file for {key}, found {matches}")
    return matches[0]


cell_tables: list[pd.DataFrame] = []
library_rows: list[dict[str, object]] = []

for gex in sorted(DATA.glob("*_GEX*_filtered_feature_bc_matrix.h5")):
    key = library_key(gex)
    adata = sc.read_10x_h5(gex)
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=[20, 50], inplace=True
    )

    cell = adata.obs[
        [
            "total_counts",
            "n_genes_by_counts",
            "pct_counts_mt",
            "pct_counts_ribo",
            "pct_counts_in_top_20_genes",
            "pct_counts_in_top_50_genes",
        ]
    ].copy()
    cell.insert(0, "barcode", cell.index.str.removesuffix("-1"))
    cell.insert(0, "library", key)
    cell_tables.append(cell.reset_index(drop=True))

    gex_barcodes = set(cell["barcode"])
    modality_stats: dict[str, int] = {}
    for modality in ("ADT", "HTO"):
        matrix_path = matching_csv(gex, modality)
        with gzip.open(matrix_path, "rt", newline="") as handle:
            header = next(csv.reader(handle))[1:]
        modality_stats[f"{modality.lower()}_all_barcodes"] = len(header)
        modality_stats[f"{modality.lower()}_matched_gex"] = len(gex_barcodes.intersection(header))

    library_rows.append(
        {
            "library": key,
            "gex_file": gex.name,
            "n_cells_filtered_matrix": adata.n_obs,
            "n_features": adata.n_vars,
            "median_umis": float(np.median(adata.obs["total_counts"])),
            "median_genes": float(np.median(adata.obs["n_genes_by_counts"])),
            "median_pct_mt": float(np.median(adata.obs["pct_counts_mt"])),
            "cells_pct_mt_gt_20": int((adata.obs["pct_counts_mt"] > 20).sum()),
            **modality_stats,
        }
    )
    print(f"Audited {key}: {adata.n_obs:,} cells")

pd.concat(cell_tables, ignore_index=True).to_csv(
    OUT / "cell_qc_unfiltered.tsv.gz", sep="\t", index=False
)
pd.DataFrame(library_rows).to_csv(OUT / "library_audit.tsv", sep="\t", index=False)
print(f"Wrote audit for {sum(r['n_cells_filtered_matrix'] for r in library_rows):,} cells")
