#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
R=Path(__file__).resolve().parents[2];O=R/'results/phase2_validation';O.mkdir(exist_ok=True)
C=['NK → Spp1 → S1pr1 → CD8','T_cell → Cd28 → Cd86 → CD8','T_cell → Lamc1 → Itga2_Itgb1 → CD4']
pd.DataFrame([{'candidate':c,'phase2a_priority':'DO_NOT_PRIORITIZE_YET','decision_reason':'Phase 2A requires direct mouse-level validation and does not alter frozen Phase 1'} for c in C]).to_csv(O/'phase2a_decision.tsv',sep='\t',index=False)
for n in ['top3_validation_summary.tsv','top3_mouse_features.tsv','top3_mouse_level_associations.tsv','top3_adjusted_models.tsv','top3_condition_robustness.tsv','top3_leave_one_mouse_out.tsv','top3_confounding_checks.tsv','top3_biological_orientation.tsv']:
 pd.DataFrame({'candidate':C,'status':['NOT ASSESSABLE']*3}).to_csv(O/n,sep='\t',index=False)
(O/'PHASE2A_TOP3_REPORT.md').write_text('# Phase 2A TOP3 computational validation\n\nPhase 1 remains frozen. Direct candidate-specific validation outputs are placeholders pending validated mouse-level feature construction; no causal priority is assigned.\n')
#!/usr/bin/env python3
raise SystemExit('Deprecated: use run_phase2a_statistical_validation.py and finalize_phase2a.py')
