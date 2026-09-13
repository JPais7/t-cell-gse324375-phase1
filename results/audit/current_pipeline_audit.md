# Current pipeline audit

## Audited state

- Starting commit: `df1200a` on `main`.
- Ranking entry point: `scripts/20_final_candidate_ranking.py`.
- Biological interaction audit: `scripts/23_interaction_biological_audit.py`.
- Mouse-removal robustness and null: `scripts/22_true_robustness.py`.
- Final integrity checks: `scripts/21_phase1_quality_control.py`.
- Principal inputs: refined T-cell object, mouse-level design/QC tables, data-driven gene evidence, TF contrasts, mechanism table with temporal support, and the external evidence registry.

## Inconsistencies found

1. Earlier LOOCV used a common mouse universe and could report more runs than candidate-informative mice.
2. Interaction-resource membership was being interpreted too readily as a biologically valid ligand–receptor direction.
3. PKM→CD44, ADAM10/17 substrate-unspecified pairs, and ITGAV in ligand position could enter the mechanistic ranking.
4. Ranking sensitivity covered only five evidence dimensions and omitted score changes.
5. Entity type and mechanistic-chain step status were not present in the evidence matrix.
6. GSE289772 was too strongly labelled `SUPPORTED`; it is contextual pharmacological support for intrinsic PKM2 biology, not gene-specific or extracellular PKM→CD44 validation.

## Corrections and regenerated outputs

- LOOCV now removes each informative mouse for each candidate, recomputes mouse-dependent evidence, score and rank, and records every fold.
- The within-candidate mouse-effect sign-permutation null uses 1,000 permutations and seed 17. It preserves effect magnitudes and candidate membership while destroying consistent direction; it is not a global pathway null.
- Interaction candidates receive explicit entity classes, validity, source, confidence and reinterpretation. Invalid/unknown pairs cannot pass the mechanistic evidence gate.
- PKM remains visible in the audit and may be considered only as an intrinsic metabolic-state hypothesis.
- Sensitivity now removes each of 13 dimensions in turn and reports rank and score changes.
- Mechanistic chain steps use `SUPPORTED`, `PARTIALLY_SUPPORTED`, `MISSING`, or `NOT_APPLICABLE`; the coherence score maps these to 1, 0.5, 0, and exclusion from the denominator.

## Provenance limitation

The manifest records the checked-out code commit at execution time. Generated result files are committed immediately after QC, so their containing results commit necessarily succeeds the code commit recorded in the manifest. No TSV was hand-edited to fabricate an analytical result; the external evidence registry is a curated input.
