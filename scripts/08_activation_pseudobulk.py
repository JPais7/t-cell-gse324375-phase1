#!/usr/bin/env python3
"""Animal-level activation contrasts without treating cells as replicates."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr, wilcoxon
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "tcells" / "GSE324375.Tcells.refined_v1.h5ad"
OUT = ROOT / "results" / "activation"
OUT.mkdir(parents=True, exist_ok=True)

RNA_SELECTOR_GENES = {
    "Cd69", "Nr4a1", "Nr4a2", "Fos", "Fosb", "Jun", "Junb", "Ier2", "Il2ra", "Tnfrsf9"
}


def bh(pvalues: np.ndarray) -> np.ndarray:
    valid = np.isfinite(pvalues)
    adjusted = np.full(len(pvalues), np.nan)
    adjusted[valid] = multipletests(pvalues[valid], method="fdr_bh")[1]
    return adjusted


def aggregate_counts(counts: sparse.spmatrix, groups: pd.Series) -> tuple[sparse.csr_matrix, pd.Index]:
    categories = pd.Index(pd.unique(groups))
    codes = pd.Categorical(groups, categories=categories).codes
    membership = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.int8), (codes, np.arange(len(codes)))),
        shape=(len(categories), len(codes)),
    )
    return (membership @ counts).tocsr(), categories


adata = ad.read_h5ad(SOURCE)
counts = sparse.csr_matrix(adata.layers["counts"])
all_results: list[pd.DataFrame] = []
pb_objects: list[ad.AnnData] = []
validation_rows: list[dict[str, object]] = []

selectors = {
    "RNA_activation": "score_activation_immediate",
    "ADT_activation": "score_activation_ADT",
}

for lineage in ("CD4", "CD8"):
    lineage_mask = adata.obs["t_lineage_provisional"].eq(lineage).to_numpy()
    lineage_obs = adata.obs.loc[lineage_mask].copy()
    lineage_counts = counts[lineage_mask]

    # Mouse-level score summaries provide an independent replication unit.
    mouse_scores = lineage_obs.groupby("mouse_id", observed=True).agg(
        n_cells=("barcode", "size"),
        RNA_activation=("score_activation_immediate", "median"),
        ADT_activation=("score_activation_ADT", "median"),
        treatment=("treatment", "first"),
        time_hours=("time_hours", "first"),
        checkpoint_blockade=("checkpoint_blockade", "first"),
    )
    mouse_scores = mouse_scores[mouse_scores["n_cells"] >= 40]
    mouse_scores.to_csv(OUT / f"{lineage}.mouse_activation_scores.tsv", sep="\t")
    rho, pvalue = spearmanr(mouse_scores["RNA_activation"], mouse_scores["ADT_activation"])
    validation_rows.append(
        {"lineage": lineage, "comparison": "mouse_RNA_vs_ADT", "n_mice": len(mouse_scores),
         "spearman_rho": rho, "pvalue": pvalue}
    )

    for selector_name, selector_column in selectors.items():
        obs = lineage_obs.copy()
        # Define high/low within each animal, preventing treatment composition from
        # determining activation status. Middle 50% is retained in the source but
        # not used in this specific contrast.
        quantiles = obs.groupby("mouse_id", observed=True)[selector_column].transform(
            lambda x: pd.Series(x).rank(pct=True).to_numpy()
        )
        obs["activation_group"] = np.select(
            [quantiles <= 0.25, quantiles >= 0.75], ["low", "high"], default="middle"
        )
        selected = obs["activation_group"].isin(["low", "high"]).to_numpy()
        selected_obs = obs.loc[selected].copy()
        selected_counts = lineage_counts[selected]
        selected_obs["pb_id"] = (
            selected_obs["mouse_id"].astype(str) + "|" + selected_obs["activation_group"]
        )
        pb_counts, pb_ids = aggregate_counts(selected_counts, selected_obs["pb_id"])
        pb_meta = (
            selected_obs.groupby("pb_id", observed=True)
            .agg(
                mouse_id=("mouse_id", "first"), activation_group=("activation_group", "first"),
                cells=("barcode", "size"), treatment=("treatment", "first"),
                time_hours=("time_hours", "first"), checkpoint_blockade=("checkpoint_blockade", "first"),
            )
            .loc[pb_ids]
        )
        keep_pb = pb_meta["cells"] >= 10
        pb_counts = pb_counts[keep_pb.to_numpy()]
        pb_meta = pb_meta.loc[keep_pb]

        paired_mice = sorted(
            set(pb_meta.loc[pb_meta.activation_group == "high", "mouse_id"])
            & set(pb_meta.loc[pb_meta.activation_group == "low", "mouse_id"])
        )
        high_index = [pb_meta.index.get_loc(f"{mouse}|high") for mouse in paired_mice]
        low_index = [pb_meta.index.get_loc(f"{mouse}|low") for mouse in paired_mice]
        library_sizes = np.asarray(pb_counts.sum(axis=1)).ravel()
        logcpm = np.log2(pb_counts.toarray() / library_sizes[:, None] * 1e6 + 0.5)
        differences = logcpm[high_index] - logcpm[low_index]
        effect = np.median(differences, axis=0)
        consistency = (differences > 0).mean(axis=0)
        stat, pvalues = wilcoxon(
            logcpm[high_index], logcpm[low_index], axis=0, zero_method="zsplit", method="approx"
        )
        result = pd.DataFrame(
            {
                "gene": adata.var_names,
                "lineage": lineage,
                "selector": selector_name,
                "n_paired_mice": len(paired_mice),
                "median_paired_log2FC": effect,
                "direction_consistency": consistency,
                "pvalue": pvalues,
                "FDR": bh(np.asarray(pvalues)),
            }
        )
        result["excluded_selector_gene"] = (
            result["gene"].isin(RNA_SELECTOR_GENES) if selector_name == "RNA_activation" else False
        )
        result["eligible_discovery_result"] = ~result["excluded_selector_gene"]
        result.sort_values(["FDR", "median_paired_log2FC"], ascending=[True, False], inplace=True)
        result.to_csv(OUT / f"{lineage}.{selector_name}.paired_DE.tsv.gz", sep="\t", index=False)
        all_results.append(result)

        pb = ad.AnnData(X=pb_counts, obs=pb_meta.copy(), var=adata.var.copy())
        pb.obs["lineage"] = lineage
        pb.obs["selector"] = selector_name
        pb_objects.append(pb)

# Cross-modal replication: a gene must have the same direction when activation is
# selected independently from RNA and surface proteins.
combined = pd.concat(all_results, ignore_index=True)
replicated_rows = []
for lineage in ("CD4", "CD8"):
    rna = combined[(combined.lineage == lineage) & (combined.selector == "RNA_activation")]
    adt = combined[(combined.lineage == lineage) & (combined.selector == "ADT_activation")]
    merged = rna.merge(adt, on=["gene", "lineage"], suffixes=("_RNA", "_ADT"))
    merged["same_direction"] = np.sign(merged["median_paired_log2FC_RNA"]) == np.sign(
        merged["median_paired_log2FC_ADT"]
    )
    merged["cross_modal_score"] = (
        -np.log10(merged[["FDR_RNA", "FDR_ADT"]].clip(lower=1e-300)).mean(axis=1)
        * np.minimum(
            merged["median_paired_log2FC_RNA"].abs(), merged["median_paired_log2FC_ADT"].abs()
        )
        * merged["same_direction"]
    )
    replicated_rows.append(merged)
pd.concat(replicated_rows, ignore_index=True).sort_values(
    "cross_modal_score", ascending=False
).to_csv(OUT / "cross_modal_activation_gene_ranking.tsv.gz", sep="\t", index=False)

ad.concat(pb_objects, join="inner", merge="same", index_unique=None).write_h5ad(
    OUT / "activation_strata_pseudobulk_counts.h5ad", compression="gzip"
)
pd.DataFrame(validation_rows).to_csv(OUT / "activation_score_mouse_validation.tsv", sep="\t", index=False)
print(pd.DataFrame(validation_rows).to_string(index=False))
