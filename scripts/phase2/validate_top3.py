#!/usr/bin/env python3
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr
ROOT=Path(__file__).resolve().parents[2]; OUT=ROOT/'results/phase2_validation'; OUT.mkdir(exist_ok=True)
C=['NK → Spp1 → S1pr1 → CD8','T_cell → Cd28 → Cd86 → CD8','T_cell → Lamc1 → Itga2_Itgb1 → CD4']
e=pd.read_csv(ROOT/'results/phase2/candidate_mouse_effects.tsv',sep='\t'); e=e[e.candidate.isin(C)].copy()
rows=[]
for c,g in e.groupby('candidate'):
 v=pd.to_numeric(g.effect_value,errors='coerce').dropna(); rho,p=(spearmanr(v,np.arange(len(v))) if len(v)>2 else (np.nan,np.nan)); rows.append({'candidate':c,'n_mice':len(v),'effect_estimate':v.mean() if len(v) else np.nan,'spearman_rho':rho,'p_value':p,'direction':'positive' if len(v) and v.mean()>0 else 'negative' if len(v) else 'NOT ASSESSABLE','primary_endpoint':'target consensus activation','association_status':'ASSOCIATION_ONLY' if len(v)>=8 else 'NOT ASSESSABLE'})
summary=pd.DataFrame(rows); summary.to_csv(OUT/'top3_validation_summary.tsv',sep='\t',index=False); e.to_csv(OUT/'top3_mouse_level_associations.tsv',sep='\t',index=False)
for name in ['top3_adjusted_models.tsv','top3_condition_robustness.tsv','top3_confounding_checks.tsv']:
 pd.DataFrame([{'candidate':c,'status':'NOT ESTIMABLE','reason':'Candidate-specific covariates unavailable; no unstable model fitted'} for c in C]).to_csv(OUT/name,sep='\t',index=False)
los=[]
for c,g in e.groupby('candidate'):
 vals=g.set_index('mouse_id').effect_value.astype(float)
 for m in vals.index: los.append({'candidate':c,'removed_mouse':m,'n_mice_remaining':len(vals)-1,'effect_mean_remaining':vals.drop(m).mean()})
pd.DataFrame(los).to_csv(OUT/'top3_leave_one_mouse_out.tsv',sep='\t',index=False)
pd.DataFrame([{'candidate':c,'phase2a_priority':'DO_NOT_PRIORITIZE_YET','decision_reason':'Targeted validation is associative and does not establish causal priority.'} for c in C]).to_csv(OUT/'phase2a_decision.tsv',sep='\t',index=False)
lines=['# Phase 2A TOP3 computational validation','','Phase 1 frozen outputs were not modified. The mouse is the inferential unit. These analyses are observational and do not establish causality.','','## Candidates']+[f'- {c}' for c in C]+['','## Decisions']+[f'- {c}: DO_NOT_PRIORITIZE_YET' for c in C]
(OUT/'PHASE2A_TOP3_REPORT.md').write_text('\n'.join(lines)+'\n')
print(f'Validated {len(C)} frozen candidates; outputs in {OUT}')
