#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad
from scipy import sparse
R=Path(__file__).resolve().parents[2]; O=R/'results/phase2_validation'; O.mkdir(exist_ok=True)
C=[('NK → Spp1 → S1pr1 → CD8','NK','Spp1','CD8',['S1pr1']),('T_cell → Cd28 → Cd86 → CD8','T_cell','Cd28','CD8',['Cd86']),('T_cell → Lamc1 → Itga2_Itgb1 → CD4','T_cell','Lamc1','CD4',['Itga2','Itgb1'])]
atlas=ad.read_h5ad(R/'results/atlas/GSE324375.atlas_v1.h5ad',backed='r'); tc=ad.read_h5ad(R/'results/tcells/GSE324375.Tcells.refined_v1.h5ad',backed='r')
ATLAS_GENE_IDX={g:atlas.var_names.get_loc(g) for g in ['Spp1','Cd28','Cd86','Lamc1'] if g in atlas.var_names}
TC_GENE_IDX={g:tc.var_names.get_loc(g) for g in ['S1pr1','Cd28','Cd86','Itga2','Itgb1','Ifng','Tnf'] if g in tc.var_names}
VECTOR_CACHE={}
def vec(a,g,layer=None):
 if g not in a.var_names:return None
 key=('atlas' if a is atlas else 'tc',g,layer)
 if key in VECTOR_CACHE:return VECTOR_CACHE[key]
 idx=(ATLAS_GENE_IDX if a is atlas else TC_GENE_IDX).get(g)
 x=a[:,idx].layers[layer] if layer else a[:,idx].X
 out=np.asarray(x.toarray() if sparse.issparse(x) else x).ravel(); VECTOR_CACHE[key]=out; return out
def one(a,mask,g):
 v=vec(a,g); ct=vec(a,g,'counts') if 'counts' in a.layers else None
 if v is None or not mask.any():return [np.nan,np.nan,np.nan,np.nan]
 vv=v[mask]; cc=ct[mask] if ct is not None else None; sm=float(cc.sum()) if cc is not None else np.nan; total=float(np.asarray(a[mask].layers['counts'].sum())) if cc is not None else np.nan; pb=np.log1p(1e6*sm/total) if total>0 else np.nan
 return [float(vv.mean()),float((vv>0).mean()),sm,float(pb)]
def feature_status(n): return 'PASS' if n>=20 else ('LOW_CELL_COUNT' if n>=10 else 'NOT_ESTIMABLE')
def safe_zscore(s):
 s=pd.to_numeric(s,errors='coerce'); o=pd.Series(np.nan,index=s.index,dtype=float); v=s.notna()&np.isfinite(s)
 if v.sum()<6 or s.loc[v].std(ddof=0)==0:return o
 o.loc[v]=(s.loc[v]-s.loc[v].mean())/s.loc[v].std(ddof=0); return o
atlas_mouse=atlas.obs.mouse_id.astype(str).to_numpy(); atlas_type=atlas.obs.provisional_cell_type.astype(str).to_numpy(); tc_mouse=tc.obs.mouse_id.astype(str).to_numpy(); tc_lineage=tc.obs.t_lineage_provisional.astype(str).to_numpy()
atlas_masks={m:atlas_mouse==m for m in np.unique(atlas_mouse)}; tc_masks={m:tc_mouse==m for m in np.unique(tc_mouse)}
atlas_type_masks={t:{m:(atlas_masks[m] & (atlas_type==t)) for m in atlas_masks} for t in np.unique(atlas_type)}
tc_lineage_masks={t:{m:(tc_masks[m] & (tc_lineage==t)) for m in tc_masks} for t in np.unique(tc_lineage)}
meta=atlas.obs.groupby(atlas.obs.mouse_id.astype(str)).agg(condition=('treatment',lambda x:';'.join(sorted(x.dropna().astype(str).unique()))),time=('time_hours',lambda x:';'.join(sorted(x.dropna().astype(str).unique()))),checkpoint_blockade=('checkpoint_blockade',lambda x:';'.join(sorted(x.dropna().astype(str).unique()))),library_or_batch=('library',lambda x:';'.join(sorted(x.dropna().astype(str).unique()))) )
rows=[]
for cand,source,sg,target,rgs in C:
 for m in sorted(meta.index):
  am=atlas_masks[m]; sm=atlas_type_masks.get(source,{}).get(m,np.zeros(len(atlas_mouse),dtype=bool)); tm=tc_lineage_masks.get(target,{}).get(m,np.zeros(len(tc_mouse),dtype=bool)); source_stats=one(atlas,sm,sg); rec=[]
  for rg in rgs: rec.append(one(tc,tm,rg))
  row={'candidate':cand,'mouse_id':m,**meta.loc[m].to_dict(),'source_population':source,'target_population':target,'source_n_cells':int(sm.sum()),'target_n_cells':int(tm.sum()),'total_relevant_cells':int(am.sum()),'source_fraction':float(sm.sum()/am.sum()) if am.sum() else np.nan,'target_fraction':float(tm.sum()/max(1,int(tc.obs.mouse_id.astype(str).eq(m).sum()))),'source_feature_status':feature_status(int(sm.sum())),'target_feature_status':feature_status(int(tm.sum()))}
  for n,v in zip(['mean','fraction_positive','sum','pseudobulk'],source_stats):row[f'{source}_{sg}_{n}']=v
  for rg,st in zip(rgs,rec):
   for n,v in zip(['mean','fraction_positive','sum','pseudobulk'],st):row[f'{target}_{rg}_{n}']=v
  for col in ['score_activation_immediate','score_activation_ADT','score_activation_consensus','score_effector_cytokine','score_cytotoxicity','score_proliferation','score_dysfunction','score_interferon_response']:
   row[f'{target}_{col.replace("score_","")}']=float(pd.to_numeric(tc.obs.loc[tm,col],errors='coerce').mean()) if tm.any() else np.nan
  row['target_Ifng']=np.nan; row['target_Tnf']=np.nan; row['target_receptor_proxy']=np.nan; row['endpoint_status']='PASS' if tm.sum()>=20 else ('LOW_CELL_COUNT' if tm.sum()>=10 else 'NOT_ESTIMABLE'); rows.append(row)
d=pd.DataFrame(rows)
for cand,source,sg,target,rgs in C:
 q=d.candidate.eq(cand); sc=f'{source}_{sg}_pseudobulk'; rc=[f'{target}_{g}_pseudobulk' for g in rgs]; zs=safe_zscore(d.loc[q,sc]); zr=sum((safe_zscore(d.loc[q,c]) for c in rc))/len(rc); d.loc[q,'axis_source_z']=zs; d.loc[q,'axis_target_z']=zr; d.loc[q,'target_receptor_proxy']=zr; d.loc[q,'candidate_axis']=(zs+zr)/2; d.loc[q,'interaction_proxy']=zs*zr
assert not d.empty and d[['candidate','mouse_id']].duplicated().sum()==0
d.to_csv(O/'top3_mouse_features.tsv',sep='\t',index=False)
pd.DataFrame([{'input_file':'results/atlas/GSE324375.atlas_v1.h5ad','level':'cell','unit':'cell','relevant_columns':'mouse_id,provisional_cell_type,treatment,time_hours,checkpoint_blockade,library','used_for':'NK/T-cell source expression and metadata'},{'input_file':'results/tcells/GSE324375.Tcells.refined_v1.h5ad','level':'cell','unit':'cell','relevant_columns':'mouse_id,t_lineage_provisional,RNA/ADT/consensus scores','used_for':'CD4/CD8 expression and endpoints'}]).to_csv(O/'phase2a_input_schema.tsv',sep='\t',index=False)
print(f'Wrote {len(d)} real candidate×mouse rows')
