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

656. fibroblast → Adam10 → Cd44 → CD8 (extracellular; score=0.744; power=HIGH)
991. neutrophil → Pdcd1lg2 → Pdcd1 → CD4 (extracellular; score=0.744; power=HIGH)
66. fibroblast → Col1a2 → Cd44 → CD4 (extracellular; score=0.742; power=HIGH)
64971. Stat4 (intrinsic; score=0.742; power=MODERATE)
64972. Atf6 (intrinsic; score=0.742; power=MODERATE)
64973. Fli1 (intrinsic; score=0.721; power=MODERATE)
17. T_cell → Pkm → Cd44 → CD8 (multicellular; score=0.777; power=HIGH)
910. NK → Tnfsf4 → Tnfrsf4 → CD8 (multicellular; score=0.775; power=HIGH)
2. T_cell → Cd80 → Ctla4 → CD4 (multicellular; score=0.764; power=HIGH)
