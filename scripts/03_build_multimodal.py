#!/usr/bin/env python3
"""Align GEX, ADT and HTO data, call hashes/doublets, and retain all cells.

This script never drops cells. It records provisional QC and doublet flags so
that thresholds can be reviewed before producing a filtered analysis object.
"""

from __future__ import annotations

import csv
import gzip
import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.mixture import GaussianMixture


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "extracted"
OUT = ROOT / "data" / "processed"
QC = ROOT / "results" / "qc"
OUT.mkdir(parents=True, exist_ok=True)
QC.mkdir(parents=True, exist_ok=True)
RUN_SCRUBLET = os.environ.get("RUN_SCRUBLET", "1") == "1"

PLANNED_HASHES = {
    "053_1": {"A0307", "A0308", "A0309", "A0310", "A0311", "A0312"},
    "053_2": {"A0307", "A0308", "A0309", "A0310", "A0311", "A0312"},
    "053_3": {"A0307", "A0308", "A0309", "A0310", "A0311", "A0312"},
    "053_6": {"A0307", "A0308"},
    "054_A": {"A0311", "A0312"},
    "054_B": {"A0311", "A0312"},
    "054_C": {"A0311", "A0312", "A0313"},
    "054_D": {"A0311", "A0312", "A0313"},
    "054_E": {"A0311", "A0312"},
}


def two_hash_gmm(counts: np.ndarray, features: list[str], index: pd.Index) -> pd.DataFrame:
    """Demultiplex sparse two-hash libraries where HashSolo is not identifiable."""
    log_counts = np.log1p(counts.astype(np.float64))
    positive = np.zeros_like(counts, dtype=bool)
    for j in range(counts.shape[1]):
        model = GaussianMixture(n_components=2, random_state=17, n_init=5)
        model.fit(log_counts[:, [j]])
        means = model.means_.ravel()
        low, high = np.argsort(means)
        # A failed tag should not be converted into a synthetic biological group.
        if means[high] - means[low] < 1.0:
            continue
        probability_high = model.predict_proba(log_counts[:, [j]])[:, high]
        positive[:, j] = probability_high >= 0.9
    n_positive = positive.sum(axis=1)
    assignment = np.full(len(index), "Negative", dtype=object)
    assignment[n_positive >= 2] = "Doublet"
    singlet = n_positive == 1
    if singlet.any():
        assignment[singlet] = np.asarray(features)[positive[singlet].argmax(axis=1)]
    result = pd.DataFrame(index=index)
    result["most_likely_hypothesis"] = np.select(
        [n_positive == 0, n_positive == 1, n_positive >= 2], [0.0, 1.0, 2.0]
    )
    result["Classification"] = assignment
    result["negative_hypothesis_probability"] = np.nan
    result["singlet_hypothesis_probability"] = np.nan
    result["doublet_hypothesis_probability"] = np.nan
    return result


def library_key(path: Path) -> str:
    name = path.name.split("_", 1)[1]
    return re.sub(r"_GEX(?:_FL)?_filtered_feature_bc_matrix\.h5$", "", name)


def matching_csv(key: str, modality: str) -> Path:
    matches = list(DATA.glob(f"*_{key}_{modality}*.csv.gz"))
    if len(matches) != 1:
        raise RuntimeError(f"Expected one {modality} file for {key}, found {matches}")
    return matches[0]


def read_selected_csv(path: Path, wanted: list[str]) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Read only wanted barcode columns from a feature-by-barcode gzipped CSV."""
    wanted_index = {barcode: i for i, barcode in enumerate(wanted)}
    with gzip.open(path, "rt", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split(",")
        source_positions: list[int] = []
        target_positions: list[int] = []
        for source_i, barcode in enumerate(header[1:]):
            target_i = wanted_index.get(barcode)
            if target_i is not None:
                source_positions.append(source_i)
                target_positions.append(target_i)
        source_positions_np = np.asarray(source_positions, dtype=np.int64)
        target_positions_np = np.asarray(target_positions, dtype=np.int64)

        features: list[str] = []
        rows: list[np.ndarray] = []
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith('"'):
                boundary = line.find('",')
                if boundary < 0:
                    raise ValueError(f"Malformed quoted feature row in {path}")
                feature = line[1:boundary].replace('""', '"')
                values = line[boundary + 2 :]
            else:
                feature, values = line.split(",", 1)
            all_values = np.fromstring(values, sep=",", dtype=np.int32)
            selected = np.zeros(len(wanted), dtype=np.int32)
            selected[target_positions_np] = all_values[source_positions_np]
            features.append(feature)
            rows.append(selected)
    return features, np.vstack(rows).T, np.isin(np.arange(len(wanted)), target_positions_np)


def clr(matrix: np.ndarray) -> np.ndarray:
    logged = np.log1p(matrix.astype(np.float32))
    return logged - logged.mean(axis=1, keepdims=True)


summary_rows: list[dict[str, object]] = []
objects: list[ad.AnnData] = []

for gex_path in sorted(DATA.glob("*_GEX*_filtered_feature_bc_matrix.h5")):
    key = library_key(gex_path)
    print(f"Processing {key}", flush=True)
    adata = sc.read_10x_h5(gex_path)
    adata.var_names_make_unique()
    adata.obs_names = pd.Index([f"{key}:{x.removesuffix('-1')}" for x in adata.obs_names])
    barcodes = [x.split(":", 1)[1] for x in adata.obs_names]
    adata.obs["library"] = key
    adata.obs["barcode"] = barcodes

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)

    adt_features, adt_counts, adt_present = read_selected_csv(matching_csv(key, "ADT"), barcodes)
    hto_features, hto_counts, hto_present = read_selected_csv(matching_csv(key, "HTO"), barcodes)
    # Feature names include slashes, which cannot be used as HDF5 keys. Store
    # matrices as arrays and keep their ordered feature names in ``uns``.
    adata.obsm["adt_counts"] = adt_counts
    adata.obsm["adt_clr"] = clr(adt_counts)
    adata.obsm["hto_counts"] = hto_counts
    adata.obsm["hto_clr"] = clr(hto_counts)
    adata.uns["adt_features"] = np.asarray(adt_features, dtype=str)
    adata.uns["hto_features"] = np.asarray(hto_features, dtype=str)
    adata.obs["adt_barcode_present"] = adt_present
    adata.obs["hto_barcode_present"] = hto_present

    # Exclude unused/background-only hashtag rows. GEO exports the full panel in
    # several libraries even when only 1–6 hashes were used. Including near-zero
    # rows destabilizes HashSolo's background model.
    hto_totals = hto_counts.sum(axis=0)
    design_key = key.removeprefix("Masopust_").rstrip("ab")
    planned = PLANNED_HASHES.get(design_key)
    active_mask = np.asarray([x in planned for x in hto_features]) if planned else (
        hto_totals >= max(100, 0.01 * float(hto_totals.max()))
    )
    active_features = [x for x, active in zip(hto_features, active_mask) if active]
    if not active_features:
        raise RuntimeError(f"No active hashtags detected for {key}")
    adata.uns["hto_active_features"] = np.asarray(active_features, dtype=str)

    # HashSolo uses the active count columns and records negative/singlet/doublet calls.
    hash_columns = [f"hash_{x}" for x in active_features]
    adata.obs[hash_columns] = hto_counts[:, active_mask]
    complete = adata.obs["hto_barcode_present"].to_numpy()
    hash_subset = adata[complete].copy()
    if len(hash_columns) == 2:
        hash_results = two_hash_gmm(
            hto_counts[complete][:, active_mask], active_features, hash_subset.obs_names
        )
        for column in hash_results.columns:
            hash_subset.obs[column] = hash_results[column]
    else:
        sc.external.pp.hashsolo(hash_subset, hash_columns)
    for column in (
        "most_likely_hypothesis",
        "negative_hypothesis_probability",
        "singlet_hypothesis_probability",
        "doublet_hypothesis_probability",
    ):
        adata.obs[column] = np.nan
        adata.obs.loc[hash_subset.obs_names, column] = pd.to_numeric(hash_subset.obs[column])
    adata.obs["Classification"] = None
    adata.obs.loc[hash_subset.obs_names, "Classification"] = hash_subset.obs["Classification"].astype(str)
    adata.obs.rename(columns={"Classification": "hto_assignment"}, inplace=True)
    adata.obs["hto_classification"] = adata.obs["most_likely_hypothesis"].map(
        {0.0: "Negative", 1.0: "Singlet", 2.0: "Doublet"}
    )
    adata.obs["hto_assignment"] = adata.obs["hto_assignment"].str.removeprefix("hash_")
    adata.obs.drop(columns=hash_columns, inplace=True)

    # Transcriptomic doublet evidence is independent of the HTO doublet call.
    try:
        if not RUN_SCRUBLET:
            raise RuntimeError("deferred (RUN_SCRUBLET=0)")
        sc.pp.scrublet(adata, expected_doublet_rate=0.08, random_state=17)
        scrublet_ok = True
    except Exception as exc:
        print(f"Scrublet failed for {key}: {exc}", flush=True)
        adata.obs["doublet_score"] = np.nan
        adata.obs["predicted_doublet"] = False
        scrublet_ok = False

    adata.obs["qc_low_complexity_provisional"] = adata.obs["n_genes_by_counts"] < 200
    adata.obs["qc_high_mito_provisional"] = adata.obs["pct_counts_mt"] > 20
    adata.obs["qc_doublet_provisional"] = (
        adata.obs["predicted_doublet"].fillna(False).astype(bool)
        | adata.obs["hto_classification"].eq("Doublet").fillna(False)
    )
    adata.obs["qc_pass_provisional"] = ~(
        adata.obs["qc_low_complexity_provisional"]
        | adata.obs["qc_high_mito_provisional"]
        | adata.obs["qc_doublet_provisional"]
    )

    summary_rows.append(
        {
            "library": key,
            "cells": adata.n_obs,
            "adt_features": len(adt_features),
            "hto_features": len(hto_features),
            "hto_active_features": ";".join(active_features),
            "adt_missing": int((~adt_present).sum()),
            "hto_missing": int((~hto_present).sum()),
            "hto_singlet": int(adata.obs["hto_classification"].eq("Singlet").sum()),
            "hto_doublet": int(adata.obs["hto_classification"].eq("Doublet").sum()),
            "hto_negative": int(adata.obs["hto_classification"].eq("Negative").sum()),
            "scrublet_ok": scrublet_ok,
            "scrublet_doublet": int(adata.obs["predicted_doublet"].fillna(False).sum()),
            "provisional_qc_pass": int(adata.obs["qc_pass_provisional"].sum()),
        }
    )
    adata.write_h5ad(OUT / f"{key}.multimodal.unfiltered.h5ad", compression="gzip")
    objects.append(adata)

pd.DataFrame(summary_rows).to_csv(QC / "multimodal_demux_summary.tsv", sep="\t", index=False)

# Concatenation preserves raw counts and makes cell identifiers globally unique.
combined = ad.concat(objects, join="inner", merge="same", index_unique=None)
combined.write_h5ad(OUT / "GSE324375.multimodal.unfiltered.h5ad", compression="gzip")
print(f"Wrote combined unfiltered object with {combined.n_obs:,} cells", flush=True)
