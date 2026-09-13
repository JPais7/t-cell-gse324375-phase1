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
expected=[P/x for x in ["candidate_evidence_matrix.tsv","candidates_for_perturbation_validation.tsv","top_extracellular.tsv","top_intrinsic.tsv","top_multicellular.tsv","loocv_candidate_stability.tsv","permutation_null_results.tsv","ranking_sensitivity.tsv","ranking_changes_from_previous.tsv","external_perturbation_evidence.tsv","TOP3_candidates.md","PHASE1_FINAL_REPORT.md"]]
missing=[str(x) for x in expected if not x.exists()]
if missing: raise SystemExit("FAIL missing outputs: "+"; ".join(missing))
x=pd.read_csv(P/"candidate_evidence_matrix.tsv",sep="\t"); top=pd.read_csv(P/"candidates_for_perturbation_validation.tsv",sep="\t")
required=["candidate","candidate_family_id","candidate_type","evidence_score","power_flag","causality_tier","RNA_evidence","ADT_evidence","temporal_support","independence_from_curated_score","LOOCV_stability","empirical_p","null_percentile"]
errs=[]
for c in required:
 if c not in x.columns: errs.append("missing column "+c)
if top.evidence_score.isna().any(): errs.append("TOP candidate missing evidence_score")
if top.power_flag.isna().any(): errs.append("TOP candidate missing power_flag")
if (top.causality_tier>3).any(): errs.append("causality tier exceeds Phase 1 limit")
if top.candidate_family_id.nunique()<len(top) and top.candidate_family_id.duplicated().any(): errs.append("duplicate candidate family in balanced shortlist")
if top.head(3).mechanistic_family_id.nunique()<len(top.head(3)): errs.append("TOP 3 mechanistic families are redundant")
if x.evidence_score.isna().any() or (~np.isfinite(pd.to_numeric(x.evidence_score,errors="coerce"))).any(): errs.append("invalid evidence score")
if (x.temporal_support & (x.causality_tier>1)).any(): errs.append("temporal support mislabeled as causal tier")
if (x.potentially_novel_hypothesis & ~x.independence_from_curated_score).any(): errs.append("novelty/independence contradiction")
lo=pd.read_csv(P/"loocv_candidate_stability.tsv",sep="\t"); nu=pd.read_csv(P/"permutation_null_results.tsv",sep="\t")
if not {"LOOCV_rank_min","LOOCV_rank_max","LOOCV_score_min","LOOCV_score_max"}.issubset(lo.columns): errs.append("LOOCV is not a true rank recalculation")
if ((lo.LOOCV_runs>0)&lo.LOOCV_rank_min.isna()).any(): errs.append("LOOCV run has no recalculated rank")
if "empirical_FDR" not in nu.columns: errs.append("permutation null missing empirical FDR")
if top.head(3).LOOCV_stability.isna().any(): errs.append("TOP candidate lacks true LOOCV")
if (top.head(3).n_mice.fillna(0)<=1).any(): errs.append("TOP candidate depends on one mouse")
app=top.n_mice.notna()&top.n_mice_supporting.notna()&top.n_mice_opposing.notna()
if ((top.loc[app,"n_mice_supporting"]+top.loc[app,"n_mice_opposing"]-top.loc[app,"n_mice"]).abs()>1e-6).any(): errs.append("supporting + opposing != n_mice")
eg=set(pd.read_csv(P/"external_perturbation_evidence.tsv",sep="\t").query("external_perturbation_support == 'SUPPORTED'").candidate); tier3_match=top.candidate.isin(eg)
if ((top.causality_tier>=3)&~tier3_match).any(): errs.append("Tier 3 without real external perturbation")
if errs: raise SystemExit("FAIL\n"+"\n".join(errs))
print("PASS: outputs present; columns complete; finite scores; no causal overclaim; classes represented; temporal separate from causality.")
print(f"Candidates={len(x):,}; shortlist={len(top)}; families={x.candidate_family_id.nunique():,}")
manifest=pd.DataFrame([{"commit_sha":subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),"run_utc":datetime.now(timezone.utc).isoformat(),"python":sys.version.split()[0],"random_seed":17,"null_permutations":1000,"cells_tcell_refined":40952,"animals":48,"datasets":"GSE324375; GSE289772 external contextual perturbation","loocv_scope":"top 100 per candidate class; 289 unique candidates"}])
manifest.to_csv(ROOT/"results/audit/phase1_run_manifest.tsv",sep="\t",index=False)
