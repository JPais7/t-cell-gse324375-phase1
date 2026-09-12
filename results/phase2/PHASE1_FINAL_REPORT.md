# Phase 1 final report

## Executive summary
This package identifies associative, animal-aware hypotheses for intratumoral T-cell activation. It preserves separate curated and unsupervised states, RNA/ADT evidence, extracellular/intrinsic/multicellular classes, power flags and explicit causal tests.

## Candidate ranking
The final score is decomposable (replication, effect, statistical strength, modality/state layers, temporal support and independence from curated score); no 80/20 LR formula or FDR-only ranking is used. All candidates are ASSOCIATION_ONLY.

## Limitations
The dataset is cross-sectional across time, lacks spatial contact and perturbation, and has underpowered relevant-vs-innate strata.

## Phase 2
Use candidates_for_perturbation_validation.tsv and TOP3_candidates.md. Test necessity and sufficiency with non-targeting, receptor/source specificity and rescue controls.

## TOP 10 candidates

1. fibroblast → Adam10 → Cd44 → CD8 (extracellular; score=0.744; power=HIGH)
2. fibroblast → Adam10 → Tspan5 → CD8 (extracellular; score=0.744; power=HIGH)
3. fibroblast → Adam10 → Tspan14 → CD8 (extracellular; score=0.744; power=HIGH)
4. Stat4 (intrinsic; score=0.492; power=MODERATE)
5. Atf6 (intrinsic; score=0.492; power=MODERATE)
6. Fli1 (intrinsic; score=0.471; power=MODERATE)
7. T_cell → Pkm → Cd44 → CD8 (multicellular; score=0.777; power=HIGH)
8. NK → Tnfsf4 → Tnfrsf4 → CD8 (multicellular; score=0.775; power=HIGH)
9. T_cell → Cd80 → Ctla4 → CD4 (multicellular; score=0.764; power=HIGH)
