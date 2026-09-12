#!/usr/bin/env python3
"""Build a discovery atlas and cluster-marker tables from analysis_v1 cells."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import harmonypy


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "processed" / "GSE324375.multimodal.analysis_v1.h5ad"
OUT = ROOT / "results" / "atlas"
OUT.mkdir(parents=True, exist_ok=True)

MARKERS = {
    "T_cell": ["Cd3d", "Cd3e", "Trbc1", "Trbc2"],
    "NK": ["Ncr1", "Klrk1", "Tyrobp", "Eomes"],
    "B_cell": ["Cd79a", "Ms4a1", "Cd37", "Cd74"],
    "monocyte_macrophage": ["Lyz2", "Csf1r", "Adgre1", "Fcgr1"],
    "dendritic": ["Flt3", "Zbtb46", "Clec9a", "Xcr1"],
    "neutrophil": ["S100a8", "S100a9", "Ly6g", "Csf3r"],
    "melanoma": ["Pmel", "Mlana", "Tyr", "Dct"],
    "endothelial": ["Pecam1", "Kdr", "Emcn", "Cdh5"],
    "fibroblast": ["Col1a1", "Col1a2", "Dcn", "Pdgfra"],
}

sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=120, frameon=False)
adata = ad.read_h5ad(SOURCE)
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, n_top_genes=4000, flavor="seurat", batch_key="library", subset=False
)

for label, genes in MARKERS.items():
    present = [gene for gene in genes if gene in adata.var_names]
    sc.tl.score_genes(adata, present, score_name=f"score_{label}", random_state=17)

sc.tl.pca(adata, n_comps=50, use_highly_variable=True, random_state=17)
sc.pp.neighbors(adata, use_rep="X_pca", key_added="unintegrated", random_state=17)
sc.tl.umap(adata, neighbors_key="unintegrated", random_state=17)
adata.obsm["X_umap_unintegrated"] = adata.obsm["X_umap"].copy()

# Harmony is used for visualization/clustering only. Counts and pseudobulk tests
# remain uncorrected and use animal-level experimental models.
harmony = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["library"], random_state=17)
# harmonypy 2.x returns cells × PCs; older Scanpy wrappers assume the previous
# orientation and transpose it incorrectly.
adata.obsm["X_pca_harmony"] = harmony.Z_corr
sc.pp.neighbors(adata, use_rep="X_pca_harmony", key_added="harmony", random_state=17)
sc.tl.umap(adata, neighbors_key="harmony", random_state=17)
sc.tl.leiden(adata, resolution=1.0, neighbors_key="harmony", key_added="leiden_harmony", random_state=17)

score_columns = [f"score_{label}" for label in MARKERS]
cluster_scores = adata.obs.groupby("leiden_harmony", observed=True)[score_columns].median()
cluster_scores["provisional_cell_type"] = cluster_scores.idxmax(axis=1).str.removeprefix("score_")
cluster_scores["manual_review_required"] = True
cluster_scores.to_csv(OUT / "cluster_marker_scores.tsv", sep="\t")
adata.obs["provisional_cell_type"] = adata.obs["leiden_harmony"].map(
    cluster_scores["provisional_cell_type"]
).astype("category")

sc.tl.rank_genes_groups(
    adata, groupby="leiden_harmony", method="t-test_overestim_var", pts=True
)
markers = sc.get.rank_genes_groups_df(adata, group=None)
markers.to_csv(OUT / "cluster_markers_all.tsv.gz", sep="\t", index=False)

for color, filename in [
    ("leiden_harmony", "umap_clusters.png"),
    ("provisional_cell_type", "umap_provisional_cell_types.png"),
    ("library", "umap_library.png"),
    ("treatment", "umap_treatment.png"),
]:
    sc.pl.umap(adata, color=color, show=False, legend_loc="on data" if "cell_type" in color else "right margin")
    plt.savefig(OUT / filename, dpi=180, bbox_inches="tight")
    plt.close("all")

adata.write_h5ad(OUT / "GSE324375.atlas_v1.h5ad", compression="gzip")
print(f"Atlas written: {adata.n_obs:,} cells, {adata.obs['leiden_harmony'].nunique()} clusters")
