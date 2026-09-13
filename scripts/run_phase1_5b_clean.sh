#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
OUT="results/phase2"
# Only downstream Phase 1.5b products are removed. Expensive upstream objects remain.
rm -f "$OUT"/{candidate_evidence_pre_audit.tsv,candidate_mouse_effects.tsv,candidate_evidence_matrix.tsv,candidates_for_perturbation_validation.tsv,top_extracellular.tsv,top_intrinsic.tsv,top_multicellular.tsv,loocv_candidate_stability.tsv,loocv_fold_details.tsv,permutation_null_results.tsv,ranking_sensitivity.tsv,ranking_changes_from_previous.tsv,interaction_biological_audit.tsv,interaction_biological_audit.md,top20_biological_audit.tsv,phase2_readiness.tsv,TOP3_candidates.md,PHASE1_FINAL_REPORT.md}
.venv/bin/python scripts/25_external_perturbation_evidence.py
.venv/bin/python scripts/24_interaction_entity_annotation.py
.venv/bin/python scripts/20_final_candidate_ranking.py
.venv/bin/python scripts/23_interaction_biological_audit.py
.venv/bin/python scripts/22_true_robustness.py
.venv/bin/python scripts/20_final_candidate_ranking.py
.venv/bin/python scripts/23_interaction_biological_audit.py
.venv/bin/python scripts/20_final_candidate_ranking.py
.venv/bin/python scripts/21_phase1_quality_control.py
