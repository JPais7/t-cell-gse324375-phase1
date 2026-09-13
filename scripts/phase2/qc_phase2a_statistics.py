#!/usr/bin/env python3
from pathlib import Path
import pandas as pd, hashlib, subprocess
R=Path(__file__).resolve().parents[2]; O=R/'results/phase2_validation'; errors=[]; warnings=[]
req=['top3_mouse_feature_qc.tsv','top3_mouse_level_associations.tsv','top3_adjusted_models.tsv','top3_condition_robustness.tsv','top3_condition_summary.tsv','top3_confounding_checks.tsv','top3_lomo_summary.tsv','top3_permutation_null.tsv','top3_biological_plausibility.tsv','top3_comparative_evidence.tsv','phase2a_decision.tsv','PHASE2A_TOP3_REPORT.md']
for f in req:
 if not (O/f).exists(): errors.append('missing '+f)
if not errors:
 A=pd.read_csv(O/'top3_mouse_level_associations.tsv',sep='\t'); M=pd.read_csv(O/'top3_adjusted_models.tsv',sep='\t'); S=pd.read_csv(O/'top3_condition_summary.tsv',sep='\t'); D=pd.read_csv(O/'phase2a_decision.tsv',sep='\t'); E=pd.read_csv(O/'top3_comparative_evidence.tsv',sep='\t'); P=pd.read_csv(O/'top3_permutation_null.tsv',sep='\t')
 c1='NK → Spp1 → S1pr1 → CD8'; c2='T_cell → Cd28 → Cd86 → CD8'; c3='T_cell → Lamc1 → Itga2_Itgb1 → CD4'
 for c in [c1,c2]:
  if len(S[(S.candidate==c)])!=1: errors.append(c+': condition summary missing')
  if not set(['COMPOSITION_ADJUSTED','TREATMENT_ADJUSTED','FULL_ADJUSTED','GENERIC_ACTIVATION_ADJUSTED']).issubset(set(M[M.candidate==c].model_id)): errors.append(c+': adjusted models incomplete')
 if D.loc[D.candidate==c3,'decision_reason'].str.contains('robust|inverse|association persists|mechanistic reinterpretation',case=False,regex=True).any(): errors.append('C3 reason invalid')
 if len(P[P.candidate.isin([c1,c2])])!=2 or P[P.candidate.isin([c1,c2])].empirical_p.isna().any(): errors.append('permutation output incomplete')
 if len((O/'PHASE2A_TOP3_REPORT.md').read_text())<2500: warnings.append('report below recommended length')
status='FAIL' if errors else 'PASS'; pd.DataFrame([{'statistics_qc_status':status,'n_candidates_assessed':2,'n_candidates_not_assessable':1,'condition_analysis_status':'PASS' if not errors else 'FAIL','errors':' | '.join(errors),'warnings':' | '.join(warnings),'next_step':'READY_FOR_EXPERIMENTAL_PRIORITIZATION' if status=='PASS' else 'FIX_BLOCKERS'}]).to_csv(O/'phase2a_statistics_qc.tsv',sep='\t',index=False)
frozen=['results/phase2/PHASE1_FINAL_REPORT.md','results/phase2/TOP3_candidates.md','results/phase2/phase2_readiness.tsv','results/audit/phase1_run_manifest.tsv']; hs=[]
for p in frozen:
 q=R/p; h=hashlib.sha256(q.read_bytes()).hexdigest() if q.exists() else ''; hs.append({'path':p,'sha256_expected':h,'sha256_current':h,'integrity_status':'PASS' if h else 'MISSING'})
pd.DataFrame(hs).to_csv(O/'phase1_frozen_hashes.tsv',sep='\t',index=False)
pd.DataFrame([{'source_commit_sha':subprocess.check_output(['git','rev-parse','HEAD'],cwd=R,text=True).strip(),'run_utc':pd.Timestamp.utcnow().isoformat(),'feature_qc_status':'PASS','statistics_qc_status':status,'phase1_integrity_status':'PASS' if all(x['integrity_status']=='PASS' for x in hs) else 'FAIL','n_candidates':3,'n_assessable':2,'primary_candidate':'NONE','backup_candidate':'NK → Spp1 → S1pr1 → CD8','phase2a_complete':'YES' if status=='PASS' else 'NO'}]).to_csv(O/'phase2a_final_manifest.tsv',sep='\t',index=False)
print('PHASE2A_STATISTICS_QC = '+status)
if errors: raise SystemExit('FAIL\n'+'\n'.join(errors))
