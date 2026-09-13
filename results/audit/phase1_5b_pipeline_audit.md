# Phase 1.5b pipeline audit

## Starting state

The run started from `main` after the Phase 1.5 commits. Upstream biological objects were preserved; only downstream candidate-audit products were regenerated. Numbered scripts 01–25 were inspected, with scripts 20–25 forming the Phase 1.5b downstream workflow.

## Previous dependency graph and problem

Previously, script 20 both constructed candidate evidence and consumed the final biological audit, LOOCV, null and final outputs. Script 23 read `candidate_evidence_matrix.tsv`, while script 22 also read that final matrix. This created a circular dependency: final ranking outputs were inputs to biological audit/robustness, which then fed the same ranking.

## Corrected acyclic graph

`upstream biological results → candidate_evidence_pre_audit.tsv → interaction_entity_annotation.tsv → interaction_biological_audit.tsv → candidate_mouse_effects.tsv → LOOCV/permutation null → final decomposable ranking/matrix → phase2_readiness → reports → QC`.

The clean runner (`scripts/run_phase1_5b_clean.sh`) deletes downstream products and rebuilds them in this order. The biological audit no longer reads the final matrix. Robustness reads the pre-audit table and biological audit, and does not read previous LOOCV/null outputs.

## Directional correction

Animal effects are stored once in `candidate_mouse_effects.tsv` and reused for supporting/opposing counts, LOOCV and sign permutations. Supporting means the expected biological sign, not the larger of positive/negative fractions. TNFSF4→TNFRSF4 has valid molecular orientation but unresolved/discordant functional direction: RNA and ADT signs differ and the mouse-level effect is slightly negative. It remains HOLD.

## Robustness correction

LOOCV uses a predeclared 279-candidate priority universe (top 100 by pre-robustness evidence per class), explicitly recorded as a priority-universe metric. Fold scores use only pre-robustness evidence plus fold-specific mouse effects; no prior LOOCV or null metrics enter the fold score. The sign null uses the same pre-robustness score, 1,000 sign flips and seed 17.

## Final status

Clean regeneration completed through QC. There are zero GO candidates, 18,185 HOLD candidates and 46,888 REJECT candidates. The leading HOLD is `NK → Tnfsf4 → Tnfrsf4 → CD8`; no TOP3 is fabricated.

## Remaining limitations

The priority universe is not all 65,073 hypotheses; this is explicit and fixed before fold results. No spatial, longitudinal or direct system-matched perturbation exists. UniProt entity annotation establishes molecular identity/localization, not physical interaction in this tumor system. Residual abundance, stress, treatment and batch confounding remain possible.
