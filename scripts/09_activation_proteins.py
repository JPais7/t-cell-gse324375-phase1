#!/usr/bin/env python3
"""Paired animal-level ADT contrasts for RNA- and ADT-defined activation."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "tcells" / "GSE324375.Tcells.refined_v1.h5ad"
OUT = ROOT / "results" / "activation"
ADT_SELECTOR = {"anti-mouse CD69", "anti-mouse CD25", "anti-mouse CD137", "anti-mouse CD134 (OX-40)", "anti-mouse CD71"}

adata = ad.read_h5ad(SOURCE)
features = [str(x) for x in adata.uns["adt_features"]]
adt = np.asarray(adata.obsm["adt_clr"], dtype=np.float32)
rows = []

for lineage in ("CD4", "CD8"):
    mask = adata.obs["t_lineage_provisional"].eq(lineage).to_numpy()
    obs = adata.obs.loc[mask].copy()
    matrix = adt[mask]
    for selector, column in {
        "RNA_activation": "score_activation_immediate",
        "ADT_activation": "score_activation_ADT",
    }.items():
        percentile = obs.groupby("mouse_id", observed=True)[column].transform(
            lambda x: pd.Series(x).rank(pct=True).to_numpy()
        )
        obs["group"] = np.select([percentile <= 0.25, percentile >= 0.75], ["low", "high"], default="middle")
        means = {}
        for (mouse, group), indices in obs.groupby(["mouse_id", "group"], observed=True).indices.items():
            if group in {"low", "high"} and len(indices) >= 10:
                means[(str(mouse), group)] = matrix[indices].mean(axis=0)
        paired = sorted({m for m, g in means if g == "high"} & {m for m, g in means if g == "low"})
        high = np.vstack([means[(m, "high")] for m in paired])
        low = np.vstack([means[(m, "low")] for m in paired])
        difference = high - low
        _, pvalue = wilcoxon(high, low, axis=0, zero_method="zsplit", method="approx")
        fdr = multipletests(pvalue, method="fdr_bh")[1]
        for i, feature in enumerate(features):
            rows.append({
                "protein": feature,
                "lineage": lineage,
                "selector": selector,
                "n_paired_mice": len(paired),
                "median_paired_CLR_difference": float(np.median(difference[:, i])),
                "direction_consistency": float((difference[:, i] > 0).mean()),
                "pvalue": float(pvalue[i]),
                "FDR": float(fdr[i]),
                "excluded_selector_protein": selector == "ADT_activation" and feature in ADT_SELECTOR,
            })

result = pd.DataFrame(rows)
result["eligible_discovery_result"] = ~result["excluded_selector_protein"]
result.sort_values(["lineage", "selector", "FDR"], inplace=True)
result.to_csv(OUT / "paired_ADT_activation_results.tsv", sep="\t", index=False)
print(result[(result.FDR < 0.05) & result.eligible_discovery_result].groupby(["lineage", "selector"]).size())
