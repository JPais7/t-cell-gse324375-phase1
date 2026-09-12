#!/usr/bin/env python3
"""Mouse-level, design-aware experimental contrasts for T-cell RNA and programs."""

from __future__ import annotations

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind, norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/tcells/GSE324375.Tcells.refined_v1.h5ad"
OUT = ROOT / "results/experimental_contrasts"
OUT.mkdir(parents=True, exist_ok=True)
MIN_CELLS = 20

PROGRAMS = [
    "score_activation_immediate", "score_activation_ADT", "score_effector_cytokine",
    "score_cytotoxicity", "score_proliferation", "score_naive_memory",
    "score_dysfunction", "score_interferon_response",
]

# Each contrast is evaluated separately inside the listed matching variables.
CONTRASTS = {
    "relevant_vs_irrelevant": ("treatment", "relevant_peptide", "irrelevant_peptide", ["time_hours", "checkpoint_blockade"]),
    "48h_vs_12h": ("time_hours", 48, 12, ["treatment", "checkpoint_blockade"]),
    "checkpoint_vs_no_checkpoint": ("checkpoint_blockade", True, False, ["treatment", "time_hours"]),
    "relevant_vs_innate": ("treatment", "relevant_peptide", "innate_agonist", ["time_hours", "checkpoint_blockade"]),
}


def bh(p):
    p = np.asarray(p, float); out = np.full(p.size, np.nan); ok = np.isfinite(p)
    if ok.any(): out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out


def combine_strata(rows, n_features):
    """Fixed-direction Stouffer meta-analysis; strata remain the replication blocks."""
    effects, ps, weights, labels, ns_a, ns_b = [], [], [], [], [], []
    for label, a, b in rows:
        effect = np.nanmean(a, axis=0) - np.nanmean(b, axis=0)
        p = ttest_ind(a, b, axis=0, equal_var=False, nan_policy="omit").pvalue
        weight = np.sqrt((a.shape[0] * b.shape[0]) / (a.shape[0] + b.shape[0]))
        effects.append(effect); ps.append(np.clip(p, 1e-300, 1)); weights.append(weight)
        labels.append(label); ns_a.append(a.shape[0]); ns_b.append(b.shape[0])
    if not rows:
        return None
    effects = np.vstack(effects); ps = np.vstack(ps); weights = np.asarray(weights)
    signs = np.sign(effects)
    z = norm.isf(ps / 2) * signs
    zmeta = np.nansum(z * weights[:, None], axis=0) / np.sqrt(np.sum(weights ** 2))
    pmeta = 2 * norm.sf(np.abs(zmeta))
    emeta = np.average(effects, axis=0, weights=weights)
    consistency = np.maximum((effects > 0).mean(0), (effects < 0).mean(0))
    return emeta, pmeta, consistency, ";".join(labels), sum(ns_a), sum(ns_b), len(rows)


adata = ad.read_h5ad(SOURCE)
keep = adata.obs.t_lineage_provisional.isin(["CD4", "CD8"]).to_numpy()
obs = adata.obs.loc[keep].copy()
counts = sparse.csr_matrix(adata.layers["counts"])[keep]

obs["pb_id"] = obs.mouse_id.astype(str) + "|" + obs.t_lineage_provisional.astype(str)
groups = pd.Index(pd.unique(obs.pb_id)); codes = pd.Categorical(obs.pb_id, categories=groups).codes
M = sparse.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(groups), len(codes)))
pb_counts = (M @ counts).tocsr(); cell_n = np.asarray(M.sum(1)).ravel()
meta = obs.groupby("pb_id", observed=True).agg(
    mouse_id=("mouse_id", "first"), lineage=("t_lineage_provisional", "first"),
    treatment=("treatment", "first"), time_hours=("time_hours", "first"),
    checkpoint_blockade=("checkpoint_blockade", "first"), n_cells=("barcode", "size"),
    **{x: (x, "median") for x in PROGRAMS},
).loc[groups]
valid = cell_n >= MIN_CELLS; meta = meta.loc[valid].copy(); pb_counts = pb_counts[valid]
lib = np.asarray(pb_counts.sum(1)).ravel()
logcpm = np.log2(pb_counts.toarray() / lib[:, None] * 1e6 + 0.5)

gene_outputs, program_outputs, stratum_gene_outputs, audit = [], [], [], []
for lineage in ["CD4", "CD8"]:
    lm = meta.lineage.eq(lineage).to_numpy()
    for cname, (factor, level_a, level_b, match) in CONTRASTS.items():
        gene_strata, program_strata = [], []
        strata = meta.loc[lm, match].drop_duplicates()
        for values in strata.itertuples(index=False, name=None):
            mask = lm.copy()
            for col, value in zip(match, values): mask &= meta[col].eq(value).to_numpy()
            ia = np.flatnonzero(mask & meta[factor].eq(level_a).to_numpy())
            ib = np.flatnonzero(mask & meta[factor].eq(level_b).to_numpy())
            label = ",".join(f"{c}={v}" for c, v in zip(match, values))
            audit.append({"lineage": lineage, "contrast": cname, "stratum": label,
                          "level_a": level_a, "level_b": level_b, "n_a": len(ia), "n_b": len(ib),
                          "included": len(ia) >= 2 and len(ib) >= 2})
            if len(ia) >= 2 and len(ib) >= 2:
                gene_strata.append((label, logcpm[ia], logcpm[ib]))
                program_strata.append((label, meta.iloc[ia][PROGRAMS].to_numpy(float), meta.iloc[ib][PROGRAMS].to_numpy(float)))
                seffect = logcpm[ia].mean(0) - logcpm[ib].mean(0)
                sp = ttest_ind(logcpm[ia], logcpm[ib], axis=0, equal_var=False).pvalue
                sf = pd.DataFrame({"feature": adata.var_names, "lineage": lineage, "contrast": cname,
                                   "level_a": level_a, "level_b": level_b, "effect_a_minus_b": seffect,
                                   "pvalue": sp, "FDR": bh(sp), "n_a": len(ia), "n_b": len(ib)})
                for col, value in zip(match, values): sf[col] = value
                stratum_gene_outputs.append(sf)
        for feature_type, rows, names, sink in [
            ("RNA", gene_strata, adata.var_names, gene_outputs),
            ("program", program_strata, PROGRAMS, program_outputs),
        ]:
            result = combine_strata(rows, len(names))
            if result is None: continue
            effect, p, consistency, labels, na, nb, nstrata = result
            frame = pd.DataFrame({"feature": names, "lineage": lineage, "contrast": cname,
                                  "level_a": level_a, "level_b": level_b, "effect_a_minus_b": effect,
                                  "pvalue": p, "FDR": bh(p), "direction_consistency": consistency,
                                  "n_a_total": na, "n_b_total": nb, "n_strata": nstrata,
                                  "included_strata": labels})
            frame["design_strength"] = np.where((nstrata >= 2) & (min(na, nb) >= 6), "supported", "fragile")
            sink.append(frame)

genes = pd.concat(gene_outputs, ignore_index=True).sort_values(["contrast", "lineage", "FDR"])
programs = pd.concat(program_outputs, ignore_index=True).sort_values(["contrast", "lineage", "FDR"])
genes.to_csv(OUT / "mouse_level_RNA_contrasts.tsv.gz", sep="\t", index=False)
programs.to_csv(OUT / "mouse_level_program_contrasts.tsv", sep="\t", index=False)
pd.DataFrame(audit).to_csv(OUT / "contrast_design_audit.tsv", sep="\t", index=False)
pd.concat(stratum_gene_outputs, ignore_index=True).to_csv(
    OUT / "stratum_RNA_contrasts.tsv.gz", sep="\t", index=False)
genes[(genes.FDR < .05) & (genes.design_strength == "supported")].to_csv(
    OUT / "supported_RNA_contrast_hits.tsv.gz", sep="\t", index=False)
print(programs.to_string(index=False, max_rows=80))
