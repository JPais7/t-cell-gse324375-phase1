#!/usr/bin/env python3
"""Infer transcription-factor activity from mouse-level T-cell pseudobulks."""

from pathlib import Path
import anndata as ad
import decoupler as dc
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import ttest_ind, norm
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/tcells/GSE324375.Tcells.refined_v1.h5ad"
OUT = ROOT / "results/regulators"; OUT.mkdir(parents=True, exist_ok=True)
MIN_CELLS = 20
CONTRASTS = {
 "relevant_vs_irrelevant": ("treatment", "relevant_peptide", "irrelevant_peptide", ["time_hours", "checkpoint_blockade"]),
 "48h_vs_12h": ("time_hours", 48, 12, ["treatment", "checkpoint_blockade"]),
 "checkpoint_vs_no_checkpoint": ("checkpoint_blockade", True, False, ["treatment", "time_hours"]),
 "relevant_vs_innate": ("treatment", "relevant_peptide", "innate_agonist", ["time_hours", "checkpoint_blockade"]),
}

def bh(x):
    x=np.asarray(x,float); o=np.full(len(x),np.nan); k=np.isfinite(x)
    if k.any(): o[k]=multipletests(x[k],method="fdr_bh")[1]
    return o

a=ad.read_h5ad(SOURCE); keep=a.obs.t_lineage_provisional.isin(["CD4","CD8"]).to_numpy()
obs=a.obs.loc[keep].copy(); cnt=sparse.csr_matrix(a.layers["counts"])[keep]
obs["pb_id"]=obs.mouse_id.astype(str)+"|"+obs.t_lineage_provisional.astype(str)
groups=pd.Index(pd.unique(obs.pb_id)); codes=pd.Categorical(obs.pb_id,categories=groups).codes
M=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(groups),len(codes)))
pc=(M@cnt).tocsr(); meta=obs.groupby("pb_id",observed=True).agg(
 mouse_id=("mouse_id","first"),lineage=("t_lineage_provisional","first"),treatment=("treatment","first"),
 time_hours=("time_hours","first"),checkpoint_blockade=("checkpoint_blockade","first"),n_cells=("barcode","size")).loc[groups]
k=meta.n_cells.to_numpy()>=MIN_CELLS; meta=meta.loc[k]; pc=pc[k]
lib=np.asarray(pc.sum(1)).ravel(); expr=np.log2(pc.toarray()/lib[:,None]*1e6+.5)
expr=pd.DataFrame(expr,index=meta.index,columns=a.var_names)
net=dc.op.dorothea(organism="mouse",levels=["A","B","C"])
net.to_csv(OUT/"dorothea_mouse_ABC_regulons.tsv.gz",sep="\t",index=False)
activity=dc.mt.ulm(expr,net,tmin=5)
if isinstance(activity,tuple): activity=activity[0]
activity.to_csv(OUT/"mouse_TF_activity.tsv.gz",sep="\t")

rows=[]
for lineage in ["CD4","CD8"]:
 lm=meta.lineage.eq(lineage).to_numpy()
 for cname,(factor,la,lb,match) in CONTRASTS.items():
  strata=[]
  for vals in meta.loc[lm,match].drop_duplicates().itertuples(index=False,name=None):
   mask=lm.copy()
   for c,v in zip(match,vals): mask &= meta[c].eq(v).to_numpy()
   ia=np.flatnonzero(mask & meta[factor].eq(la).to_numpy()); ib=np.flatnonzero(mask & meta[factor].eq(lb).to_numpy())
   if len(ia)>=2 and len(ib)>=2:
    ef=activity.iloc[ia].mean(0).to_numpy()-activity.iloc[ib].mean(0).to_numpy()
    p=ttest_ind(activity.iloc[ia],activity.iloc[ib],axis=0,equal_var=False).pvalue
    w=np.sqrt(len(ia)*len(ib)/(len(ia)+len(ib))); strata.append((ef,p,w,len(ia),len(ib)))
  if not strata: continue
  ef=np.vstack([x[0] for x in strata]); p=np.vstack([x[1] for x in strata]); w=np.array([x[2] for x in strata])
  z=norm.isf(np.clip(p,1e-300,1)/2)*np.sign(ef); zm=np.nansum(z*w[:,None],0)/np.sqrt((w*w).sum())
  pm=2*norm.sf(abs(zm)); em=np.average(ef,axis=0,weights=w)
  f=pd.DataFrame({"TF":activity.columns,"lineage":lineage,"contrast":cname,"level_a":la,"level_b":lb,
                  "activity_effect_a_minus_b":em,"pvalue":pm,"FDR":bh(pm),"n_strata":len(strata),
                  "direction_consistency":np.maximum((ef>0).mean(0),(ef<0).mean(0))})
  rows.append(f)
res=pd.concat(rows,ignore_index=True).sort_values(["contrast","lineage","FDR"])
res.to_csv(OUT/"TF_activity_experimental_contrasts.tsv",sep="\t",index=False)
print(res.groupby(["contrast","lineage"]).head(5).to_string(index=False))
