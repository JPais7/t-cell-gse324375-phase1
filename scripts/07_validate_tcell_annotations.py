#!/usr/bin/env python3
"""Generate RNA–ADT concordance and T-cell annotation QC tables."""

from pathlib import Path

import anndata as ad
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "tcells" / "GSE324375.Tcells.refined_v1.h5ad"
DOUBLETS = ROOT / "results" / "tcells" / "Tcells.with_doublet_calls.h5ad"
OUT = ROOT / "results" / "tcells"

adata = ad.read_h5ad(SOURCE, backed="r")
obs = adata.obs.copy()

pairs = {
    "CD4_RNA_vs_ADT": ("score_lineage_CD4", "ADT_CD4"),
    "CD8_RNA_vs_ADT_a": ("score_lineage_CD8", "ADT_CD8a"),
    "CD8_RNA_vs_ADT_b": ("score_lineage_CD8", "ADT_CD8b"),
    "activation_RNA_vs_ADT": ("score_activation_immediate", "score_activation_ADT"),
    "dysfunction_RNA_vs_PD1": ("score_dysfunction", "ADT_PD1"),
    "dysfunction_RNA_vs_TIM3": ("score_dysfunction", "ADT_TIM3"),
}

rows = []
for comparison, (left, right) in pairs.items():
    rho, pvalue = spearmanr(obs[left], obs[right], nan_policy="omit")
    rows.append({"comparison": comparison, "spearman_rho": rho, "pvalue": pvalue, "n_cells": len(obs)})
pd.DataFrame(rows).to_csv(OUT / "rna_adt_concordance.tsv", sep="\t", index=False)

pd.crosstab(obs["tcell_cluster"], obs["t_lineage_provisional"]).to_csv(
    OUT / "cluster_by_lineage.tsv", sep="\t"
)
pd.crosstab(obs["t_lineage_provisional"], obs["congenic_provisional"]).to_csv(
    OUT / "lineage_by_congenic.tsv", sep="\t"
)
obs.groupby(["treatment", "time_hours", "checkpoint_blockade"], observed=True).agg(
    cells=("barcode", "size"),
    activation_RNA_median=("score_activation_immediate", "median"),
    activation_ADT_median=("score_activation_ADT", "median"),
    cytotoxicity_median=("score_cytotoxicity", "median"),
    dysfunction_median=("score_dysfunction", "median"),
).reset_index().to_csv(OUT / "programs_by_condition_descriptive.tsv", sep="\t", index=False)

dbl = ad.read_h5ad(DOUBLETS, backed="r").obs
dbl.groupby("library", observed=True).agg(
    cells=("barcode", "size"),
    predicted_doublets=("predicted_doublet", "sum"),
    median_doublet_score=("doublet_score", "median"),
).to_csv(OUT / "scrublet_by_library.tsv", sep="\t")

print(pd.DataFrame(rows).to_string(index=False))
