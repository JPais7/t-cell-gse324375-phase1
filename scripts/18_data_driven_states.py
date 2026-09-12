#!/usr/bin/env python3
"""Unsupervised T-cell states built without curated activation genes."""

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
ad.settings.allow_write_nullable_strings = True

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/tcells/GSE324375.Tcells.refined_v1.h5ad"
OUT=ROOT/"results/data_driven"; OUT.mkdir(parents=True,exist_ok=True)

CURATED={"Cd69","Nr4a1","Nr4a2","Fos","Fosb","Jun","Junb","Ier2","Il2ra","Tnfrsf9",
         "Ifng","Tnf","Ccl3","Ccl4","Ccl5","Xcl1","Xcl2","Nkg7","Prf1","Gzma","Gzmb","Gzmk","Ctsw",
         "Mki67","Top2a","Stmn1","Tyms","Ube2c","Cenpf","Pdcd1","Havcr2","Lag3","Tox","Tigit","Entpd1","Cd160",
         "Stat1","Isg15","Ifit1","Ifit2","Ifit3","Irf7","Usp18"}

a=ad.read_h5ad(SOURCE)
a=a[a.obs.t_lineage_provisional.isin(["CD4","CD8"])].copy()
X=sparse.csr_matrix(a.layers["counts"])
a.X=X
sc.pp.normalize_total(a,target_sum=1e4); sc.pp.log1p(a)
a.var["exclude_curated"]=[g in CURATED for g in a.var_names]
a.var["exclude_qc"]=(a.var.mt.fillna(False)|a.var.ribo.fillna(False))
a.var["feature_for_state"]=(~a.var.exclude_curated)&(~a.var.exclude_qc)
sub=a[:,a.var.feature_for_state].copy()
sc.pp.highly_variable_genes(sub,n_top_genes=2500,flavor="cell_ranger",batch_key="library",subset=True)
sc.pp.scale(sub,max_value=10)
sc.tl.pca(sub,n_comps=40,svd_solver="arpack",random_state=17)
sc.pp.neighbors(sub,n_neighbors=20,n_pcs=30,random_state=17)
sc.tl.leiden(sub,resolution=0.6,key_added="data_driven_cluster",random_state=17,flavor="igraph",directed=False)
sc.tl.umap(sub,random_state=17)

# Identify the activation-enriched unsupervised state only after clustering, using
# curated scores as an external interpretation (never as clustering features).
sub.obs["curated_activation_external"]=a.obs.loc[sub.obs_names,"score_activation_immediate"].to_numpy()
state_means=sub.obs.groupby("data_driven_cluster",observed=True).agg(
 n_cells=("barcode","size"), mean_curated=("curated_activation_external","mean"),
 median_curated=("curated_activation_external","median"),
 mean_ADT=("score_activation_ADT","mean"),
).sort_values("mean_curated",ascending=False)
state_means["data_driven_activation_state"]=False
if len(state_means): state_means.iloc[0,state_means.columns.get_loc("data_driven_activation_state")]=True
state_means.to_csv(OUT/"state_summary.tsv",sep="\t")

sub.obs["data_driven_activation_state"] = sub.obs.data_driven_cluster.eq(state_means.index[0])
sub.obs["data_driven_activation_state"] = sub.obs.data_driven_activation_state.astype(bool)
sub.obs[["data_driven_cluster","data_driven_activation_state","curated_activation_external"]].to_csv(OUT/"cell_state_assignments.tsv.gz",sep="\t")
sub.write_h5ad(OUT/"Tcell_data_driven_state_atlas.h5ad",compression="gzip")

# Markers of the selected state, with all curated genes removed from discovery.
sc.tl.rank_genes_groups(sub,"data_driven_cluster",groups=[state_means.index[0]],reference="rest",method="wilcoxon",pts=True)
markers=sc.get.rank_genes_groups_df(sub,group=state_means.index[0])
markers=markers[~markers.names.isin(CURATED)].copy()
markers.rename(columns={"names":"gene","logfoldchanges":"log2FC_approx","pvals_adj":"FDR"},inplace=True)
markers.to_csv(OUT/"data_driven_state_markers.tsv.gz",sep="\t",index=False)

# Explicit circularity/audit table requested in the specification.
cross=pd.read_csv(ROOT/"results/activation/cross_modal_activation_gene_ranking.tsv.gz",sep="\t")
cross=cross.groupby(["gene","lineage"],observed=True).agg(cross_modal_score=("cross_modal_score","max"),
  RNA_FDR=("FDR_RNA","min"),ADT_FDR=("FDR_ADT","min")).reset_index()
state=markers[["gene","FDR","log2FC_approx"]].drop_duplicates("gene").rename(columns={"FDR":"data_driven_FDR","log2FC_approx":"data_driven_log2FC"})
audit=cross.merge(state,on="gene",how="outer")
audit["used_in_curated_score"]=audit.gene.isin(CURATED)
audit["associated_with_data_driven_state"]=audit.data_driven_FDR<0.05
audit["independent_candidate"]=(~audit.used_in_curated_score)&audit.associated_with_data_driven_state
audit.to_csv(OUT/"gene_evidence_circularity_audit.tsv.gz",sep="\t",index=False)
print(state_means.to_string())
print(f"Data-driven markers: {(markers.FDR<0.05).sum():,}; independent candidates: {audit.independent_candidate.sum():,}")
