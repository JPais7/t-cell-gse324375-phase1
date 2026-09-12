#!/usr/bin/env python3
"""Attach animal/condition metadata and create a reversible QC view."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "processed" / "GSE324375.multimodal.unfiltered.h5ad"
OUT = ROOT / "data" / "processed"
QC = ROOT / "results" / "qc"


# Mapping reconstructed from the GEO feature README. Keys are library suffixes;
# values map active HTO identifiers to (treatment, time_hours, checkpoint_blockade).
DESIGN: dict[str, dict[str, tuple[str, int, bool]]] = {
    "053_1": {
        **{f"A030{x}": ("baseline", 0, False) for x in (7, 8, 9)},
        **{f"A03{x}": ("relevant_peptide", 12, False) for x in (10, 11, 12)},
    },
    "053_2": {
        **{f"A030{x}": ("relevant_peptide", 48, False) for x in (7, 8, 9)},
        **{f"A03{x}": ("irrelevant_peptide", 12, False) for x in (10, 11, 12)},
    },
    "053_3": {
        **{f"A030{x}": ("irrelevant_peptide", 48, False) for x in (7, 8, 9)},
        **{f"A03{x}": ("relevant_peptide", 12, True) for x in (10, 11, 12)},
    },
    "053_6": {f"A030{x}": ("innate_agonist", 48, True) for x in (7, 8)},
    "054_A": {f"A03{x}": ("relevant_peptide", 12, True) for x in (11, 12)},
    "054_B": {f"A03{x}": ("relevant_peptide", 48, True) for x in (11, 12)},
    "054_C": {f"A03{x}": ("irrelevant_peptide", 12, True) for x in (11, 12, 13)},
    "054_D": {f"A03{x}": ("irrelevant_peptide", 48, True) for x in (11, 12, 13)},
    "054_E": {f"A03{x}": ("innate_agonist", 12, True) for x in (11, 12)},
}


def design_key(library: str) -> str:
    return library.removeprefix("Masopust_").rstrip("ab")


adata = ad.read_h5ad(SOURCE)
treatments: list[str | None] = []
times: list[int | None] = []
icb: list[bool | None] = []
mouse_ids: list[str | None] = []

for library, assignment, classification in zip(
    adata.obs["library"], adata.obs["hto_assignment"], adata.obs["hto_classification"]
):
    if classification != "Singlet" or pd.isna(assignment):
        treatments.append(None); times.append(None); icb.append(None); mouse_ids.append(None)
        continue
    details = DESIGN.get(design_key(str(library)), {}).get(str(assignment))
    if details is None:
        treatments.append(None); times.append(None); icb.append(None); mouse_ids.append(None)
        continue
    treatment, time_hours, checkpoint = details
    treatments.append(treatment); times.append(time_hours); icb.append(checkpoint)
    mouse_ids.append(f"{library}:{assignment}")

adata.obs["treatment"] = pd.Categorical(treatments)
adata.obs["time_hours"] = pd.array(times, dtype="Int64")
adata.obs["checkpoint_blockade"] = pd.array(icb, dtype="boolean")
adata.obs["mouse_id"] = pd.Categorical(mouse_ids)
adata.obs["design_resolved"] = adata.obs["mouse_id"].notna()

# This is an explicit, reviewable view—not deletion from the source object.
adata.obs["analysis_eligible_v1"] = (
    adata.obs["design_resolved"]
    & adata.obs["hto_classification"].eq("Singlet")
    & (adata.obs["n_genes_by_counts"] >= 200)
    & (adata.obs["pct_counts_mt"] <= 20)
)

adata.obs.to_csv(QC / "cell_manifest_with_design.tsv.gz", sep="\t")
(
    adata.obs.groupby(
        ["treatment", "time_hours", "checkpoint_blockade", "mouse_id"],
        observed=True,
    )
    .agg(cells=("barcode", "size"), eligible=("analysis_eligible_v1", "sum"))
    .reset_index()
    .to_csv(QC / "experimental_design_cell_counts.tsv", sep="\t", index=False)
)

adata.write_h5ad(SOURCE, compression="gzip")
eligible = adata[adata.obs["analysis_eligible_v1"]].copy()
eligible.write_h5ad(OUT / "GSE324375.multimodal.analysis_v1.h5ad", compression="gzip")
print(f"Design resolved: {adata.obs['design_resolved'].sum():,}/{adata.n_obs:,}")
print(f"Analysis v1: {eligible.n_obs:,}/{adata.n_obs:,} cells")
