#!/usr/bin/env python3
"""T-cell-specific doublet detection, RNA/ADT validation and state atlas."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import harmonypy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import zscore


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "atlas" / "GSE324375.atlas_v1.h5ad"
REFERENCE = ROOT / "data" / "processed" / "Masopust_053_1a.multimodal.unfiltered.h5ad"
OUT = ROOT / "results" / "tcells"
OUT.mkdir(parents=True, exist_ok=True)

PROGRAMS = {
    "activation_immediate": ["Cd69", "Nr4a1", "Nr4a2", "Fos", "Fosb", "Jun", "Junb", "Ier2", "Il2ra", "Tnfrsf9"],
    "effector_cytokine": ["Ifng", "Tnf", "Ccl3", "Ccl4", "Ccl5", "Xcl1", "Xcl2"],
    "cytotoxicity": ["Nkg7", "Prf1", "Gzma", "Gzmb", "Gzmk", "Ctsw"],
    "proliferation": ["Mki67", "Top2a", "Stmn1", "Tyms", "Ube2c", "Cenpf"],
    "naive_memory": ["Tcf7", "Ccr7", "Lef1", "Sell", "Il7r", "Mal"],
    "dysfunction": ["Pdcd1", "Havcr2", "Lag3", "Tox", "Tigit", "Entpd1", "Cd160"],
    "interferon_response": ["Stat1", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Irf7", "Usp18"],
}

LINEAGE_RNA = {
    "CD4": ["Cd4", "Il7r", "Ltb", "Mal"],
    "CD8": ["Cd8a", "Cd8b1", "Ccl5", "Nkg7"],
}

ADT_NAMES = {
    "CD4": "anti-mouse CD4",
    "CD8a": "anti-mouse CD8a",
    "CD8b": "anti-mouse CD8b (Ly-3)",
    "CD3": "anti-mouse CD3",
    "CD45_1": "anti-mouse CD45.1",
    "CD45_2": "anti-mouse CD45.2",
    "CD69": "anti-mouse CD69",
    "CD25": "anti-mouse CD25",
    "CD137": "anti-mouse CD137",
    "OX40": "anti-mouse CD134 (OX-40)",
    "PD1": "anti-mouse CD279 (PD-1)",
    "TIM3": "anti-mouse CD366 (Tim-3)",
    "LAG3": "anti-mouse CD223 (LAG-3)",
    "KLRG1": "anti-mouse/human KLRG1 (MAFA)",
    "CD44": "anti-mouse/human CD44",
    "CD62L": "anti-mouse CD62L",
    "CD71": "anti-mouse CD71",
}


def safe_z(values: pd.Series | np.ndarray) -> np.ndarray:
    result = zscore(np.asarray(values, dtype=float), nan_policy="omit")
    return np.nan_to_num(result, nan=0.0, posinf=0.0, neginf=0.0)


atlas = ad.read_h5ad(SOURCE)
reference = ad.read_h5ad(REFERENCE, backed="r")
adt_features = [str(x) for x in reference.uns["adt_features"]]
feature_index = {feature: i for i, feature in enumerate(adt_features)}

tcells = atlas[atlas.obs["provisional_cell_type"].eq("T_cell")].copy()
tcells.uns["adt_features"] = np.asarray(adt_features, dtype=str)

# Add interpretable ADT columns without altering the full 129-feature matrix.
for short, full in ADT_NAMES.items():
    if full in feature_index:
        tcells.obs[f"ADT_{short}"] = tcells.obsm["adt_clr"][:, feature_index[full]]

# Detect same-hash/transcriptomic doublets within each library. HTO doublets have
# already been excluded from analysis_v1.
tcells.X = tcells.layers["counts"].copy()
sc.pp.scrublet(tcells, batch_key="library", expected_doublet_rate=0.08, random_state=17)
tcells.obs["tcell_analysis_eligible"] = ~tcells.obs["predicted_doublet"].fillna(False)
tcells.write_h5ad(OUT / "Tcells.with_doublet_calls.h5ad", compression="gzip")

t = tcells[tcells.obs["tcell_analysis_eligible"]].copy()
sc.pp.normalize_total(t, target_sum=1e4)
sc.pp.log1p(t)
sc.pp.highly_variable_genes(t, n_top_genes=3000, flavor="seurat", batch_key="library")

for label, genes in {**PROGRAMS, **{f"lineage_{k}": v for k, v in LINEAGE_RNA.items()}}.items():
    present = [gene for gene in genes if gene in t.var_names]
    sc.tl.score_genes(t, present, score_name=f"score_{label}", random_state=17)

# RNA/ADT lineage calls: agreement is preferred, disagreement remains ambiguous.
cd4_evidence = safe_z(t.obs["score_lineage_CD4"]) + safe_z(t.obs["ADT_CD4"])
cd8_adt = (t.obs["ADT_CD8a"] + t.obs["ADT_CD8b"]) / 2
cd8_evidence = safe_z(t.obs["score_lineage_CD8"]) + safe_z(cd8_adt)
difference = cd4_evidence - cd8_evidence
t.obs["t_lineage_provisional"] = pd.Categorical(
    np.select([difference >= 0.75, difference <= -0.75], ["CD4", "CD8"], default="ambiguous")
)

# Congenic phenotype is reported without assuming which compartment is OT-I.
congenic_difference = safe_z(t.obs["ADT_CD45_1"]) - safe_z(t.obs["ADT_CD45_2"])
t.obs["congenic_provisional"] = pd.Categorical(
    np.select(
        [congenic_difference >= 1.0, congenic_difference <= -1.0],
        ["CD45.1_high", "CD45.2_high"],
        default="ambiguous",
    )
)

adt_activation = np.column_stack(
    [safe_z(t.obs[f"ADT_{name}"]) for name in ("CD69", "CD25", "CD137", "OX40", "CD71")]
).mean(axis=1)
t.obs["score_activation_ADT"] = adt_activation
t.obs["score_activation_consensus"] = (
    safe_z(t.obs["score_activation_immediate"]) + safe_z(adt_activation)
) / 2

sc.tl.pca(t, n_comps=50, mask_var="highly_variable", random_state=17)
harmony = harmonypy.run_harmony(t.obsm["X_pca"], t.obs, ["library"], random_state=17)
t.obsm["X_pca_harmony"] = harmony.Z_corr
sc.pp.neighbors(t, use_rep="X_pca_harmony", random_state=17)
sc.tl.umap(t, random_state=17)
sc.tl.leiden(
    t,
    resolution=1.0,
    key_added="tcell_cluster",
    flavor="igraph",
    n_iterations=2,
    directed=False,
    random_state=17,
)

program_columns = [f"score_{name}" for name in PROGRAMS]
cluster_programs = t.obs.groupby("tcell_cluster", observed=True)[program_columns].median()
cluster_programs["dominant_program_provisional"] = (
    cluster_programs.idxmax(axis=1).str.removeprefix("score_")
)
cluster_programs["manual_review_required"] = True
cluster_programs.to_csv(OUT / "tcell_cluster_programs.tsv", sep="\t")
t.obs["dominant_program_provisional"] = t.obs["tcell_cluster"].map(
    cluster_programs["dominant_program_provisional"]
).astype("category")

sc.tl.rank_genes_groups(t, groupby="tcell_cluster", method="t-test_overestim_var", pts=True)
sc.get.rank_genes_groups_df(t, group=None).to_csv(
    OUT / "tcell_cluster_markers.tsv.gz", sep="\t", index=False
)

summary = (
    t.obs.groupby(
        ["tcell_cluster", "t_lineage_provisional", "dominant_program_provisional"], observed=True
    )
    .size()
    .rename("cells")
    .reset_index()
)
summary.to_csv(OUT / "tcell_cluster_summary.tsv", sep="\t", index=False)

for color, filename in [
    ("tcell_cluster", "umap_tcell_clusters.png"),
    ("t_lineage_provisional", "umap_tcell_lineage.png"),
    ("dominant_program_provisional", "umap_tcell_programs.png"),
    ("score_activation_consensus", "umap_activation_consensus.png"),
    ("score_dysfunction", "umap_dysfunction.png"),
    ("treatment", "umap_tcell_treatment.png"),
]:
    sc.pl.umap(t, color=color, show=False)
    plt.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close("all")

t.write_h5ad(OUT / "GSE324375.Tcells.refined_v1.h5ad", compression="gzip")
print(f"T cells before Scrublet: {tcells.n_obs:,}")
print(f"Predicted transcriptomic doublets: {tcells.obs['predicted_doublet'].sum():,}")
print(f"Refined T-cell atlas: {t.n_obs:,} cells, {t.obs['tcell_cluster'].nunique()} clusters")
