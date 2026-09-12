#!/usr/bin/env python3
"""Sensitivity analysis of activation-associated gene effects under alternate definitions."""

from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/"results/tcells/GSE324375.Tcells.refined_v1.h5ad"
BASE=ROOT/"results/activation/cross_modal_activation_gene_ranking.tsv.gz"
MECH=ROOT/"results/mechanisms/mechanism_shortlist_20.tsv"
OUT=ROOT/"results/sensitivity"; OUT.mkdir(parents=True,exist_ok=True)

a=ad.read_h5ad(SOURCE); counts=sparse.csr_matrix(a.layers["counts"])
base=pd.read_csv(BASE,sep="\t")
genes=set(base.sort_values("cross_modal_score",ascending=False).groupby("lineage",observed=True).head(250).gene)
if MECH.exists():
 m=pd.read_csv(MECH,sep="\t"); genes.update(m.ligand); genes.update(g for x in m.receptor for g in str(x).split("_"))
genes=sorted(genes & set(a.var_names)); gi=a.var_names.get_indexer(genes)

variants=[
 ("primary_q25",.25,20,False),("q20",.20,20,False),("q33",.33,20,False),
 ("min40",.25,40,False),("min80",.25,80,False),("exclude_congenic_ambiguous",.25,20,True),
]
effects=[]
for lineage in ["CD4","CD8"]:
 for selector,col in [("RNA_activation","score_activation_immediate"),("ADT_activation","score_activation_ADT")]:
  lineage_mask=a.obs.t_lineage_provisional.eq(lineage).to_numpy()
  for vname,q,mincells,exclude_congenic in variants:
   mask=lineage_mask.copy()
   if exclude_congenic: mask &= ~a.obs.congenic_provisional.eq("ambiguous").to_numpy()
   obs=a.obs.loc[mask].copy(); X=counts[mask][:,gi]
   ranks=obs.groupby("mouse_id",observed=True)[col].transform(lambda x:x.rank(pct=True))
   obs["grp"]=np.select([ranks<=q,ranks>=1-q],["low","high"],default="mid")
   use=obs.grp.ne("mid").to_numpy(); obs=obs.loc[use]; X=X[use]
   obs["pb"]=obs.mouse_id.astype(str)+"|"+obs.grp
   groups=pd.Index(pd.unique(obs.pb)); codes=pd.Categorical(obs.pb,categories=groups).codes
   M=sparse.csr_matrix((np.ones(len(codes)),(codes,np.arange(len(codes)))),shape=(len(groups),len(codes)))
   pc=(M@X).toarray(); meta=obs.groupby("pb",observed=True).agg(mouse=("mouse_id","first"),grp=("grp","first"),cells=("barcode","size")).loc[groups]
   mice=sorted(set(meta.loc[(meta.grp=="high")&(meta.cells>=mincells),"mouse"]) & set(meta.loc[(meta.grp=="low")&(meta.cells>=mincells),"mouse"]))
   dif=[]
   for mouse in mice:
    hi=meta.index.get_loc(f"{mouse}|high"); lo=meta.index.get_loc(f"{mouse}|low")
    h=np.log2(pc[hi]/pc[hi].sum()*1e6+.5); l=np.log2(pc[lo]/pc[lo].sum()*1e6+.5); dif.append(h-l)
   if not dif: continue
   effect=np.median(np.vstack(dif),0)
   effects.append(pd.DataFrame({"gene":genes,"lineage":lineage,"selector":selector,"variant":vname,"effect":effect,"n_paired_mice":len(mice)}))
e=pd.concat(effects,ignore_index=True); e.to_csv(OUT/"activation_effects_by_variant.tsv.gz",sep="\t",index=False)
summary=[]
for (lin,sel),x in e.groupby(["lineage","selector"],observed=True):
 pivot=x.pivot(index="gene",columns="variant",values="effect")
 for v in pivot.columns:
  rho,p=spearmanr(pivot.primary_q25,pivot[v]); concord=(np.sign(pivot.primary_q25)==np.sign(pivot[v])).mean()
  summary.append({"lineage":lin,"selector":sel,"variant":v,"spearman_vs_primary":rho,"pvalue":p,"direction_concordance":concord})
pd.DataFrame(summary).to_csv(OUT/"activation_sensitivity_summary.tsv",sep="\t",index=False)
print(pd.DataFrame(summary).to_string(index=False))
