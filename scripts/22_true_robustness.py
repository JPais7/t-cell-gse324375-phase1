#!/usr/bin/env python3
"""True mouse-removal stability and stratified sign-permutation null for priority candidates."""
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/phase2"; SEED=17; NPERM=1000
ev=pd.read_csv(OUT/"candidate_evidence_matrix.tsv",sep="\t")
priority=pd.concat([pd.read_csv(OUT/f"top_{c}.tsv",sep="\t").head(100) for c in ["extracellular","intrinsic","multicellular"]],ignore_index=True).drop_duplicates("candidate")

# Recover per-mouse intrinsic high-low effects from the stored RNA-selector pseudobulks.
pb=ad.read_h5ad(ROOT/"results/activation/activation_strata_pseudobulk_counts.h5ad")
pmask=pb.obs.selector.eq("RNA_activation").to_numpy(); pm=pb.obs.loc[pmask].copy(); pc=sparse.csr_matrix(pb.X)[pmask]
lib=np.asarray(pc.sum(1)).ravel(); log=np.log2(pc.toarray()/lib[:,None]*1e6+.5)
per_candidate={}
for row in priority.itertuples(index=False):
 vals={}
 if row.candidate_type=="intrinsic" and row.candidate in pb.var_names:
  j=pb.var_names.get_loc(row.candidate)
  for mouse in pm.mouse_id.astype(str).unique():
   q=pm.mouse_id.astype(str).eq(mouse)&pm.lineage.eq(row.target)&pm.activation_group.isin(["high","low"])
   z=np.flatnonzero(q.to_numpy())
   if len(z)==2:
    hi=z[pm.iloc[z].activation_group.to_numpy()=="high"]; lo=z[pm.iloc[z].activation_group.to_numpy()=="low"]
    if len(hi)==1 and len(lo)==1: vals[mouse]=float(log[hi[0],j]-log[lo[0],j])
 per_candidate[row.candidate]=vals

# Recover source-ligand expression and target activation per mouse for LR hypotheses.
atlas=ad.read_h5ad(ROOT/"results/atlas/GSE324375.atlas_v1.h5ad")
tc=ad.read_h5ad(ROOT/"results/tcells/GSE324375.Tcells.refined_v1.h5ad",backed="r")
act={lin:tc.obs[tc.obs.t_lineage_provisional.eq(lin)].groupby("mouse_id",observed=True).score_activation_immediate.median() for lin in ["CD4","CD8"]}
for row in priority[priority.candidate_type.ne("intrinsic")].itertuples(index=False):
 vals={}
 if isinstance(row.ligand,str) and row.ligand in atlas.var_names:
  mask=atlas.obs.provisional_cell_type.eq(row.source).to_numpy(); x=sparse.csr_matrix(atlas.X)[mask,atlas.var_names.get_loc(row.ligand)].toarray().ravel(); o=atlas.obs.loc[mask,["mouse_id","treatment","time_hours","checkpoint_blockade"]].copy(); o["x"]=x
  g=o.groupby("mouse_id",observed=True).agg(x=("x","mean"),n_source_cells=("x","size"),condition=("treatment","first"),time=("time_hours","first"),icb=("checkpoint_blockade","first")); g=g[g.n_source_cells>=20]; y=act.get(row.target,pd.Series(dtype=float)); common=g.index.intersection(y.index)
  if len(common)>=4:
   xx=g.loc[common,"x"]; yy=y.loc[common]; strata=g.loc[common].condition.astype(str)+"|"+g.loc[common].time.astype(str)+"|"+g.loc[common].icb.astype(str); xr=xx-xx.groupby(strata).transform("mean"); yr=yy-yy.groupby(strata).transform("mean")
   vals={str(m):float(a*b) for m,a,b in zip(common,xr,yr)}
 per_candidate[row.candidate]=vals

def score(row,values):
 v=np.asarray(list(values.values()),float)
 if len(v)<2:return np.nan
 consistency=max((v>0).mean(),(v<0).mean()); repl=min(len(v)/10,1)*consistency
 constants=[row.effect_component,row.statistical_component,row.RNA_component,row.ADT_component,row.state_component,row.temporal_component,row.interaction_component,row.stability_component,row.null_model_component,row.biological_validity_component,row.mechanistic_evidence_coherence,row.independence_component,row.experimental_testability_component]
 value=float(np.mean([repl]+constants))
 if getattr(row,"external_perturbation_support","NOT_ASSESSED")=="SUPPORTED": value=float(np.mean([repl]+constants+[row.external_perturbation_component]))
 return value

universe=sorted(set().union(*[set(v) for v in per_candidate.values()])); rows=[]; folds=[]; rng=np.random.default_rng(SEED); nullrows=[]
base=priority.set_index("candidate")
rank_records={c:[] for c in base.index}; score_records={c:[] for c in base.index}; fail={c:0 for c in base.index}
# Each focal candidate runs exactly once per mouse informative for that candidate.
for focal in base.index:
 for mouse in per_candidate.get(focal,{}):
  scores={}
  for c,row in base.iterrows():
   vals={m:v for m,v in per_candidate.get(c,{}).items() if m!=mouse}; s=score(row,vals)
   if np.isfinite(s): scores[c]=s
  ranks=pd.Series(scores).rank(ascending=False,method="min")
  if focal in scores:
   rank_records[focal].append(float(ranks[focal])); score_records[focal].append(scores[focal]); folds.append({"candidate":focal,"removed_mouse":mouse,"fold_score":scores[focal],"fold_rank":ranks[focal],"eligible_after_removal":True,"failure_reason":""})
  else:
   fail[focal]+=1; folds.append({"candidate":focal,"removed_mouse":mouse,"fold_score":np.nan,"fold_rank":np.nan,"eligible_after_removal":False,"failure_reason":"insufficient remaining mice"})
for c,row in base.iterrows():
 vals=per_candidate.get(c,{}); rr=np.asarray(rank_records[c]); ss=np.asarray(score_records[c]); n=len(vals)
 vv=np.asarray(list(vals.values()),float); support=max(int((vv>0).sum()),int((vv<0).sum())) if len(vv) else 0
 rows.append({"candidate":c,"candidate_family_id":row.candidate_family_id,"n_mice":n,"n_mice_supporting":support,"n_mice_opposing":n-support,"direction_consistency":support/n if n else np.nan,"mouse_effect_median":float(np.median(vv)) if n else np.nan,"mouse_effect_mean":float(np.mean(vv)) if n else np.nan,"mouse_effect_IQR":float(np.subtract(*np.percentile(vv,[75,25]))) if n else np.nan,"mouse_effect_sign":"positive" if n and np.median(vv)>0 else ("negative" if n and np.median(vv)<0 else "undetermined"),"LOOCV_runs":len(rr),"LOOCV_failed_runs":fail[c],"LOOCV_top1_runs":int((rr<=1).sum()),"LOOCV_top3_runs":int((rr<=3).sum()),"LOOCV_top10_runs":int((rr<=10).sum()),"LOOCV_top1_fraction":float((rr<=1).mean()) if len(rr) else np.nan,"LOOCV_top10_fraction":float((rr<=10).mean()) if len(rr) else np.nan,"LOOCV_top3_fraction":float((rr<=3).mean()) if len(rr) else np.nan,"LOOCV_rank_median":float(np.median(rr)) if len(rr) else np.nan,"LOOCV_rank_mean":float(np.mean(rr)) if len(rr) else np.nan,"LOOCV_rank_min":float(np.min(rr)) if len(rr) else np.nan,"LOOCV_rank_max":float(np.max(rr)) if len(rr) else np.nan,"LOOCV_rank_IQR":float(np.subtract(*np.percentile(rr,[75,25]))) if len(rr) else np.nan,"LOOCV_score_median":float(np.median(ss)) if len(ss) else np.nan,"LOOCV_score_min":float(np.min(ss)) if len(ss) else np.nan,"LOOCV_score_max":float(np.max(ss)) if len(ss) else np.nan,"LOOCV_failure_reason":"insufficient remaining mice" if fail[c] else ""})
 v=np.asarray(list(vals.values()),float); observed=score(row,vals)
 if len(v)>=2:
  ns=[]
  for _ in range(NPERM): ns.append(score(row,{str(i):x*s for i,(x,s) in enumerate(zip(v,rng.choice([-1,1],len(v))))}))
  ns=np.asarray(ns); nullrows.append({"candidate":c,"candidate_family_id":row.candidate_family_id,"observed_score":observed,"null_mean":ns.mean(),"null_sd":ns.std(),"null_95":np.quantile(ns,.95),"null_99":np.quantile(ns,.99),"empirical_p":(1+(ns>=observed).sum())/(NPERM+1),"null_percentile":(ns<=observed).mean(),"n_permutations":NPERM,"null_model_type":"within-candidate mouse-effect sign permutation","random_seed":SEED})
 else: nullrows.append({"candidate":c,"candidate_family_id":row.candidate_family_id,"observed_score":observed,"null_model_type":"NOT_APPLICABLE","n_permutations":0,"random_seed":SEED})
lo=pd.DataFrame(rows); nu=pd.DataFrame(nullrows); ok=nu.empirical_p.notna(); nu["empirical_FDR"]=np.nan
if ok.any(): nu.loc[ok,"empirical_FDR"]=multipletests(nu.loc[ok,"empirical_p"],method="fdr_bh")[1]
lo.to_csv(OUT/"loocv_candidate_stability.tsv",sep="\t",index=False); pd.DataFrame(folds).to_csv(OUT/"loocv_fold_details.tsv",sep="\t",index=False); nu.rename(columns={"null_95":"null_p95","null_99":"null_p99"}).to_csv(OUT/"permutation_null_results.tsv",sep="\t",index=False)
print(f"Priority candidates={len(priority)}; mouse universe={len(universe)}; true LOOCV rows={len(lo)}; permutation rows={len(nu)}")
