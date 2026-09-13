#!/usr/bin/env python3
"""Fail-fast integrity checks for the final Phase 1 handoff."""
from pathlib import Path
import hashlib
import subprocess
import sys
from datetime import datetime, timezone
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; P=ROOT/"results/phase2"
expected=[P/x for x in ["candidate_evidence_pre_audit.tsv","candidate_mouse_effects.tsv","candidate_evidence_matrix.tsv","candidates_for_perturbation_validation.tsv","top_extracellular.tsv","top_intrinsic.tsv","top_multicellular.tsv","loocv_candidate_stability.tsv","loocv_fold_details.tsv","permutation_null_results.tsv","ranking_sensitivity.tsv","ranking_changes_from_previous.tsv","external_perturbation_evidence.tsv","interaction_entity_annotation.tsv","interaction_biological_audit.tsv","interaction_biological_audit.md","top20_biological_audit.tsv","phase2_readiness.tsv","TOP3_candidates.md","PHASE1_FINAL_REPORT.md"]]+[ROOT/"results/audit/phase1_5b_pipeline_audit.md"]
missing=[str(x) for x in expected if not x.exists()]
if missing: raise SystemExit("FAIL missing outputs: "+"; ".join(missing))
x=pd.read_csv(P/"candidate_evidence_matrix.tsv",sep="\t"); top=pd.read_csv(P/"candidates_for_perturbation_validation.tsv",sep="\t")
required=["candidate","candidate_family_id","mechanistic_family_id","candidate_type","entity_type","mechanism_level","evidence_score","power_flag","causality_tier","RNA_evidence","ADT_evidence","RNA_direction","ADT_direction","RNA_ADT_direction_concordant","expected_direction","direction_mechanistically_concordant","temporal_support","independence_from_curated_score","LOOCV_stability","empirical_p","empirical_FDR","null_percentile","biological_validity_component","mechanistic_evidence_coherence","experimental_testability_component"]
errs=[]
for c in required:
 if c not in x.columns: errs.append("missing column "+c)
if top.evidence_score.isna().any(): errs.append("TOP candidate missing evidence_score")
if top.power_flag.isna().any(): errs.append("TOP candidate missing power_flag")
if (top.causality_tier>=4).any() or (x.causality_tier>=4).any(): errs.append("Tier 4 appears in Phase 1.5")
if top.candidate_family_id.nunique()<len(top) and top.candidate_family_id.duplicated().any(): errs.append("duplicate candidate family in balanced shortlist")
if top.head(3).mechanistic_family_id.nunique()<len(top.head(3)): errs.append("TOP 3 mechanistic families are redundant")
if x.evidence_score.isna().any() or (~np.isfinite(pd.to_numeric(x.evidence_score,errors="coerce"))).any(): errs.append("invalid evidence score")
# Temporal evidence may coexist with Tier 2, but cannot itself create Tier 3.
if ((x.causality_tier>=3)&~x.external_perturbation_support.eq("SUPPORTED")).any(): errs.append("causal tier incompatible with external evidence")
if (x.potentially_novel_hypothesis & ~x.independence_from_curated_score).any(): errs.append("novelty/independence contradiction")
lo=pd.read_csv(P/"loocv_candidate_stability.tsv",sep="\t"); nu=pd.read_csv(P/"permutation_null_results.tsv",sep="\t")
effects=pd.read_csv(P/"candidate_mouse_effects.tsv",sep="\t")
if effects.empty or not {"candidate","mouse_id","effect_value","expected_direction","effect_supports_hypothesis"}.issubset(effects.columns): errs.append("missing single-source candidate mouse effects")
if not {"LOOCV_rank_min","LOOCV_rank_max","LOOCV_score_min","LOOCV_score_max"}.issubset(lo.columns): errs.append("LOOCV is not a true rank recalculation")
if ((lo.LOOCV_runs>0)&lo.LOOCV_rank_min.isna()).any(): errs.append("LOOCV run has no recalculated rank")
if ((lo.LOOCV_runs+lo.LOOCV_failed_runs)!=lo.n_mice).any(): errs.append("LOOCV folds do not match informative mice")
if "empirical_FDR" not in nu.columns: errs.append("permutation null missing empirical FDR")
if top.head(3).LOOCV_stability.isna().any(): errs.append("TOP candidate lacks true LOOCV")
if (top.head(3).n_mice.fillna(0)<=1).any(): errs.append("TOP candidate depends on one mouse")
app=top.n_mice.notna()&top.n_mice_supporting.notna()&top.n_mice_opposing.notna()
if ((top.loc[app,"n_mice_supporting"]+top.loc[app,"n_mice_opposing"]-top.loc[app,"n_mice"]).abs()>1e-6).any(): errs.append("supporting + opposing != n_mice")
eg=set(pd.read_csv(P/"external_perturbation_evidence.tsv",sep="\t").query("external_perturbation_support == 'SUPPORTED'").candidate); tier3_match=top.candidate.isin(eg)
if ((top.causality_tier>=3)&~tier3_match).any(): errs.append("Tier 3 without real external perturbation")
if "valid_for_mechanistic_ranking" in x and (~x.loc[x.candidate.isin(top.candidate),"valid_for_mechanistic_ranking"].fillna(False)).any(): errs.append("biologically invalid interaction entered TOP candidates")
bio=pd.read_csv(P/"interaction_biological_audit.tsv",sep="\t")
bad=bio.ligand_entity_type.str.contains("intracellular|metabolic|transcription factor",case=False,na=False)
if bio.loc[bad,"interaction_validity"].isin(["VALID_LR","MEMBRANE_CONTACT","SHEDDING_PROCESSING"]).any(): errs.append("intracellular/metabolic ligand retained as valid LR")
sens=pd.read_csv(P/"ranking_sensitivity.tsv",sep="\t")
needed={"replication","effect","FDR/statistics","RNA","ADT","data_driven_state","interaction","temporal","LOOCV","permutation_null","biological_validity","mechanistic_coherence","external_perturbation","independence"}
if not needed.issubset(set(sens.omitted_component)): errs.append("ranking sensitivity missing dimensions")
if not {"rank_without_component","score_change"}.issubset(sens.columns): errs.append("ranking sensitivity lacks rank/score changes")
if top[["suggested_perturbation","expected_readout","negative_control","specificity_control"]].isna().any().any(): errs.append("final candidate has no interpretable experimental test")
ready=pd.read_csv(P/"phase2_readiness.tsv",sep="\t"); go=ready.phase2_decision.eq("GO")
if ((go)&(~ready.biological_status.eq("PASS")|~ready.statistical_status.eq("PASS")|~ready.experimental_status.eq("PASS")|~ready.LOOCV_status.eq("PASS")|~ready.null_status.eq("PASS"))).any(): errs.append("Phase 2 GO despite failed critical gate")
if set(top.candidate)!=set(ready.loc[go,"candidate"]): errs.append("shortlist and Phase 2 GO decisions disagree")
for script,output in [(ROOT/"scripts/20_final_candidate_ranking.py",P/"candidate_evidence_matrix.tsv"),(ROOT/"scripts/22_true_robustness.py",P/"loocv_candidate_stability.tsv"),(ROOT/"scripts/23_interaction_biological_audit.py",P/"interaction_biological_audit.tsv"),(ROOT/"scripts/24_interaction_entity_annotation.py",P/"interaction_entity_annotation.tsv")]:
 if output.stat().st_mtime<script.stat().st_mtime: errs.append(f"stale output older than generator: {output.name}")
source23=(ROOT/"scripts/23_interaction_biological_audit.py").read_text(); source22=(ROOT/"scripts/22_true_robustness.py").read_text()
if "candidate_evidence_pre_audit.tsv" not in source23 or "candidate_evidence_matrix.tsv" in source23: errs.append("biological audit has circular final-matrix dependency")
if "candidate_evidence_pre_audit.tsv" not in source22 or "candidate_evidence_matrix.tsv" in source22: errs.append("robustness has circular final-matrix dependency")
if "pd.read_csv(OUT/\"candidate_evidence_matrix.tsv\"" in source22 or "pd.read_csv(OUT/\"loocv_candidate_stability.tsv\"" in source22 or "pd.read_csv(OUT/\"permutation_null_results.tsv\"" in source22: errs.append("robustness reads prior robustness outputs")
if not (ROOT/"scripts/run_phase1_5b_clean.sh").exists(): errs.append("missing clean-run workflow")
if errs: raise SystemExit("FAIL\n"+"\n".join(errs))
print("PASS: Phase 1.5 outputs are current, finite, biologically gated, mouse-robust, sensitivity-complete and experimentally interpretable.")
print(f"Candidates={len(x):,}; shortlist={len(top)}; families={x.candidate_family_id.nunique():,}")
mouse_universe=pd.read_csv(P/"loocv_fold_details.tsv",sep="\t").removed_mouse.nunique()
manifest=pd.DataFrame([{"commit_sha":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"run_utc":datetime.now(timezone.utc).isoformat(),"python":sys.version.split()[0],"random_seed":17,"n_candidates":len(x),"n_mice":int(mouse_universe),"n_permutations":int(nu.n_permutations.max()),"LOOCV_candidates":len(lo),"input_datasets":"GSE324375","external_datasets":"GSE289772 (contextual); GSE314342 (not assessed)","pipeline_version":"Phase1.5 UniProt biological-audit v3"}])
manifest.to_csv(ROOT/"results/audit/phase1_run_manifest.tsv",sep="\t",index=False)
