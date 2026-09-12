#!/usr/bin/env python3
"""Fail-fast integrity checks for the final Phase 1 handoff."""
from pathlib import Path
import hashlib
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; P=ROOT/"results/phase2"
expected=[P/x for x in ["candidate_evidence_matrix.tsv","candidates_for_perturbation_validation.tsv","top_extracellular.tsv","top_intrinsic.tsv","top_multicellular.tsv","TOP3_candidates.md","PHASE1_FINAL_REPORT.md"]]
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
if top.candidate_type.nunique()<3: errs.append("candidate classes not fairly represented")
if x.evidence_score.isna().any() or (~np.isfinite(pd.to_numeric(x.evidence_score,errors="coerce"))).any(): errs.append("invalid evidence score")
if (x.temporal_support & (x.causality_tier>1)).any(): errs.append("temporal support mislabeled as causal tier")
if (x.potentially_novel_hypothesis & ~x.independence_from_curated_score).any(): errs.append("novelty/independence contradiction")
if errs: raise SystemExit("FAIL\n"+"\n".join(errs))
print("PASS: outputs present; columns complete; finite scores; no causal overclaim; classes represented; temporal separate from causality.")
print(f"Candidates={len(x):,}; shortlist={len(top)}; families={x.candidate_family_id.nunique():,}")
