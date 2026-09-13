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
ev=pd.read_csv(OUT/"candidate_evidence_pre_audit.tsv",sep="\t",low_memory=False)
bio=pd.read_csv(OUT/"interaction_biological_audit.tsv",sep="\t",low_memory=False)[["candidate","interaction_validity","valid_for_mechanistic_ranking"]].drop_duplicates("candidate")
ev=ev.merge(bio,on="candidate",how="left"); ev["biological_validity_component"]=np.select([ev.candidate_type.eq("intrinsic"),ev.interaction_validity.eq("VALID_LR"),ev.interaction_validity.isin(["MEMBRANE_CONTACT","SHEDDING_PROCESSING"])],[1.0,1.0,.8],default=0.0)
ev["mechanistic_evidence_coherence"]=np.where(ev.candidate_type.eq("intrinsic"),ev.RNA_evidence.astype(float),ev.interaction_evidence.astype(float))
ev["experimental_testability_component"]=1.0
# Predeclared priority universe: top 100 by pre-robustness evidence per candidate class,
# constructed before any fold/null result exists.
ev["pre_robustness_score"]=ev[[c for c in ["effect_size","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score"] if c in ev]].apply(pd.to_numeric,errors="coerce").fillna(0).mean(axis=1)
ev["effect_component"]=(pd.to_numeric(ev.effect_size,errors="coerce").abs()/pd.to_numeric(ev.effect_size,errors="coerce").abs().quantile(.95)).clip(0,1).fillna(0); ev["statistical_component"]=(-np.log10(pd.to_numeric(ev.FDR,errors="coerce").fillna(1).clip(1e-300,1))/10).clip(0,1)
priority=pd.concat([ev[ev.candidate_type.eq(c)].sort_values("pre_robustness_score",ascending=False).head(100) for c in ["extracellular","intrinsic","multicellular"]],ignore_index=True).drop_duplicates("candidate")

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

# Single source of truth for all animal-level robustness calculations.
effect_rows=[]
for row in priority.itertuples(index=False):
 expected=getattr(row,"expected_direction","UNRESOLVED")
 for mouse,value in per_candidate.get(row.candidate,{}).items():
  support=(value>0) if expected=="POSITIVE_ACTIVATION" else ((value<0) if expected=="NEGATIVE_ACTIVATION" else False)
  effect_rows.append({"candidate":row.candidate,"candidate_type":row.candidate_type,"mouse_id":mouse,"effect_value":value,"effect_definition":"intrinsic high-minus-low RNA log2 expression" if row.candidate_type=="intrinsic" else "within-stratum centered source expression multiplied by centered target activation","expected_direction":expected,"effect_supports_hypothesis":support if expected!="UNRESOLVED" else "UNRESOLVED","condition":None,"time":None,"checkpoint_blockade":None,"n_source_cells":None,"target_lineage":row.target_lineage})
pd.DataFrame(effect_rows).to_csv(OUT/"candidate_mouse_effects.tsv",sep="\t",index=False)

def score(row,values):
 v=np.asarray(list(values.values()),float)
 if len(v)<2:return np.nan
 expected=getattr(row,"expected_direction","UNRESOLVED"); support=(v>0) if expected=="POSITIVE_ACTIVATION" else ((v<0) if expected=="NEGATIVE_ACTIVATION" else np.zeros(len(v),dtype=bool)); repl=min(len(v)/10,1)*(support.mean() if expected!="UNRESOLVED" else 0.0)
 constants=[row.effect_component,row.statistical_component,float(row.RNA_evidence),float(row.ADT_evidence),float(row.data_driven_state_evidence),float(row.temporal_support),float(row.interaction_evidence),row.biological_validity_component,row.mechanistic_evidence_coherence,float(row.independence_from_curated_score),row.experimental_testability_component]
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
   rank_records[focal].append(float(ranks[focal])); score_records[focal].append(scores[focal]); focal_vals={m:v for m,v in per_candidate.get(focal,{}).items() if m!=mouse}; fv=np.asarray(list(focal_vals.values()),float); fr=getattr(base.loc[focal],"expected_direction","UNRESOLVED"); fs=int((fv>0).sum()) if fr=="POSITIVE_ACTIVATION" else (int((fv<0).sum()) if fr=="NEGATIVE_ACTIVATION" else 0); fo=int((fv<0).sum()) if fr=="POSITIVE_ACTIVATION" else (int((fv>0).sum()) if fr=="NEGATIVE_ACTIVATION" else 0); folds.append({"candidate":focal,"removed_mouse":mouse,"n_mice_remaining":len(fv),"n_supporting_remaining":fs,"n_opposing_remaining":fo,"direction_consistency_remaining":fs/len(fv) if fr!="UNRESOLVED" and len(fv) else np.nan,"median_effect_remaining":float(np.median(fv)) if len(fv) else np.nan,"fold_pre_robustness_score":scores[focal],"fold_score":scores[focal],"fold_rank":ranks[focal],"n_competing_candidates":len(scores),"eligible_after_removal":True,"failure_reason":""})
  else:
   fail[focal]+=1; folds.append({"candidate":focal,"removed_mouse":mouse,"n_mice_remaining":len(vals)-1,"n_supporting_remaining":np.nan,"n_opposing_remaining":np.nan,"direction_consistency_remaining":np.nan,"median_effect_remaining":np.nan,"fold_pre_robustness_score":np.nan,"fold_score":np.nan,"fold_rank":np.nan,"n_competing_candidates":len(scores),"eligible_after_removal":False,"failure_reason":"insufficient remaining mice"})
for c,row in base.iterrows():
 vals=per_candidate.get(c,{}); rr=np.asarray(rank_records[c]); ss=np.asarray(score_records[c]); n=len(vals)
 vv=np.asarray(list(vals.values()),float); expected=getattr(row,"expected_direction","UNRESOLVED"); pos=int((vv>0).sum()); neg=int((vv<0).sum()); support=pos if expected=="POSITIVE_ACTIVATION" else (neg if expected=="NEGATIVE_ACTIVATION" else 0); opposing=neg if expected=="POSITIVE_ACTIVATION" else (pos if expected=="NEGATIVE_ACTIVATION" else 0)
 rows.append({"candidate":c,"candidate_family_id":row.candidate_family_id,"n_mice":n,"n_mice_supporting":support,"n_mice_opposing":opposing,"direction_consistency":support/n if expected!="UNRESOLVED" and n else np.nan,"expected_direction":expected,"n_positive":pos,"n_negative":neg,"n_zero_or_undetermined":n-pos-neg,"mouse_effect_median":float(np.median(vv)) if n else np.nan,"mouse_effect_mean":float(np.mean(vv)) if n else np.nan,"mouse_effect_IQR":float(np.subtract(*np.percentile(vv,[75,25]))) if n else np.nan,"mouse_effect_sign":"positive" if n and np.median(vv)>0 else ("negative" if n and np.median(vv)<0 else "undetermined"),"LOOCV_runs":len(rr),"LOOCV_failed_runs":fail[c],"LOOCV_top1_runs":int((rr<=1).sum()),"LOOCV_top3_runs":int((rr<=3).sum()),"LOOCV_top10_runs":int((rr<=10).sum()),"LOOCV_top1_fraction":float((rr<=1).mean()) if len(rr) else np.nan,"LOOCV_top10_fraction":float((rr<=10).mean()) if len(rr) else np.nan,"LOOCV_top3_fraction":float((rr<=3).mean()) if len(rr) else np.nan,"LOOCV_rank_median":float(np.median(rr)) if len(rr) else np.nan,"LOOCV_rank_mean":float(np.mean(rr)) if len(rr) else np.nan,"LOOCV_rank_min":float(np.min(rr)) if len(rr) else np.nan,"LOOCV_rank_max":float(np.max(rr)) if len(rr) else np.nan,"LOOCV_rank_IQR":float(np.subtract(*np.percentile(rr,[75,25]))) if len(rr) else np.nan,"LOOCV_score_median":float(np.median(ss)) if len(ss) else np.nan,"LOOCV_score_min":float(np.min(ss)) if len(ss) else np.nan,"LOOCV_score_max":float(np.max(ss)) if len(ss) else np.nan,"LOOCV_failure_reason":"insufficient remaining mice" if fail[c] else ""})
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
