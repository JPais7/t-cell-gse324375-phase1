# Phase 1.5d final candidate-aware pipeline audit

This audit records the final candidate-aware clean run. The run comprises pre-audit candidate construction, entity annotation, upstream biological audit, animal-level effects, independent LOOCV, directional/two-sided sign nulls, final ranking, post-ranking summary and QC.

The authoritative run metadata and counts are recorded in `phase1_run_manifest.tsv`. Outputs are frozen only after QC reports PASS.

## Run identity and dependency graph

Pipeline version: Phase1.5d final candidate-aware clean run; random seed 17; 1,000 sign permutations. The acyclic workflow is: upstream results → external perturbation evidence → entity annotation → candidate pre-audit → biological audit → candidate mouse effects → LOOCV and permutation null → final ranking → post-ranking summary → QC.

## Freeze status

The leading hypothesis remains a HOLD unless unchanged gates produce GO. No causal claim, therapeutic benefit, or artificial TOP3 is inferred. The manifest is the authoritative record of commit SHA, timestamp, candidate counts, null classes, LOOCV universe and QC status.
