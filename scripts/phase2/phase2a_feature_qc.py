#!/usr/bin/env python3
from pathlib import Path
import numpy as np
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
summ=[]
for c,g in d.groupby('candidate'):
 endpoint=next((col for col in g.columns if col.endswith('_activation_consensus')),None)
 cc=np.isfinite(pd.to_numeric(g.candidate_axis,errors='coerce'))&np.isfinite(pd.to_numeric(g[endpoint],errors='coerce')) if endpoint else np.zeros(len(g),dtype=bool)
 summ.append({'candidate':c,'n_total_mice':len(g),'n_source_pass':int((g.source_feature_status=='PASS').sum()),'n_source_low_count':int((g.source_feature_status=='LOW_CELL_COUNT').sum()),'n_source_not_estimable':int((g.source_feature_status=='NOT_ESTIMABLE').sum()),'n_target_pass':int((g.target_feature_status=='PASS').sum()),'n_target_low_count':int((g.target_feature_status=='LOW_CELL_COUNT').sum()),'n_target_not_estimable':int((g.target_feature_status=='NOT_ESTIMABLE').sum()),'n_endpoint_pass':int((g.endpoint_status=='PASS').sum()),'n_endpoint_low_count':int((g.endpoint_status=='LOW_CELL_COUNT').sum()),'n_endpoint_not_estimable':int((g.endpoint_status=='NOT_ESTIMABLE').sum()),'n_axis_complete':int(cc.sum()),'n_complete_cases':int(cc.sum()),'primary_predictor':'candidate_axis','primary_endpoint':'target_activation_consensus','candidate_feature_status':'READY_FOR_STATISTICAL_VALIDATION' if cc.sum()>=6 else 'NOT_ASSESSABLE'})
pd.DataFrame(summ).to_csv(O/'top3_mouse_feature_qc.tsv',sep='\t',index=False)
if errs: raise SystemExit('FAIL\n'+'\n'.join(errs))
print(f'FEATURE QC PASS: {len(d)} rows, {d.mouse_id.nunique()} mice, 3 candidates')
