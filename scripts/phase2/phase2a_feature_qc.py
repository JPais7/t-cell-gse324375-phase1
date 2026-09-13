#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2]; O=R/'results/phase2_validation'; f=O/'top3_mouse_features.tsv'; d=pd.read_csv(f,sep='\t')
required={'candidate','mouse_id','candidate_axis','interaction_proxy'}; missing=required-set(d.columns); activation_cols=[c for c in d.columns if c.endswith('_activation_consensus')];
if not activation_cols: missing.add('target activation consensus')
errs=[]
if d.empty: errs.append('feature table is empty')
if missing: errs.append('missing columns: '+','.join(sorted(missing)))
if d[['candidate','mouse_id']].duplicated().any(): errs.append('duplicate candidate×mouse rows')
if d.candidate.nunique()!=3: errs.append('expected exactly three frozen candidates')
if d.mouse_id.nunique()<8: errs.append('insufficient informative mice')
status='PASS' if not errs else 'FAIL'
pd.DataFrame([{'feature_qc_status':status,'n_rows':len(d),'n_candidates':d.candidate.nunique(),'n_mice':d.mouse_id.nunique(),'missing_required_columns':';'.join(sorted(missing)),'errors':';'.join(errs),'next_step':'READY_FOR_STATISTICAL_VALIDATION' if status=='PASS' else 'FIX_FEATURE_TABLE'}]).to_csv(O/'phase2a_feature_qc.tsv',sep='\t',index=False)
if errs: raise SystemExit('FAIL\n'+'\n'.join(errs))
print(f'FEATURE QC PASS: {len(d)} rows, {d.mouse_id.nunique()} mice, 3 candidates')
