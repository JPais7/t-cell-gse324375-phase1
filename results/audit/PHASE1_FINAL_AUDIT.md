# Phase 1 final scientific and computational audit

## Scope and pipeline

The repository is a reproducible numbered pipeline: GEO manifest and QC (`01–04`), atlas and T-cell refinement (`05–07`), activation RNA/ADT pseudobulk and pathways (`08–10`), source–ligand–receptor and LIANA evidence (`11–12`), legacy mechanism ranking (`13`), experimental contrasts/temporal checks/TF activity/sensitivity (`14–17`), unsupervised state discovery and balanced candidate classes (`18–19`). The raw archive, extracted matrices, virtual environment and `.h5ad` objects are excluded from Git; metadata, scripts, figures and tabular results are versioned.

## Already correct

- HTO-demultiplexed mouse/library-hash units are the biological replicate for mouse-level inference.
- RNA and ADT remain separate evidence dimensions; selector features are excluded from corresponding discovery contrasts.
- Raw counts and reversible QC flags are retained; transcriptomic and HTO doublet evidence are distinct.
- Comparisons are stratified by time and checkpoint status where design permits.
- LIANA and ligand–receptor scores are expression plausibility/association evidence, never causal evidence.

## Corrections made

- Added unsupervised states built without curated activation/effector/proliferation/dysfunction/interferon genes.
- Added animal-level state fractions, informative-mouse counts and power flags; cell-level state markers are labeled exploratory.
- Replaced the binary “novelty” interpretation with `independent_from_curated_score`, `known_activation_gene`, `mechanistically_supported` and `potentially_novel_hypothesis` fields in the final ranking.
- Added balanced extracellular, intrinsic and multicellular candidate classes and a decomposable evidence score; the legacy 80/20 LR score is retained only as an input feature, not the final scientific rank.
- Added explicit animal replication, effect, consistency, modality, state, temporal and power components, plus minimum evidence gates for TOP 3.
- Added Phase 2 candidate table, TOP 3 rationale and final report.
- Added candidate-family collapsing, leave-one-mouse-out stability, a 1,000-permutation rank null model, and fail-fast quality control (`scripts/21_phase1_quality_control.py`). The null model is a robustness benchmark for the composite score, not a replacement for molecular FDR.

## What remains associative

All Phase 1 candidates are associative. A positive RNA/ADT association, a source–receptor compatible pair, a data-driven state association or cross-sectional temporal ordering does not establish necessity, sufficiency, physical contact or signaling direction. TF activities are inferred regulatory scores, not direct TF measurements.

## What cannot be inferred from GSE324375 alone

The data do not provide same-animal longitudinal trajectories, spatial contact, direct protein measurements for every gene, perturbational effects, or proof that a ligand is delivered to a receptor. The relevant-vs-innate comparison is underpowered in several strata. No Phase 1 result can be called causal or therapeutic.

No external perturbational or literature dataset was supplied/integrated; `external_perturbation_support` and `literature_novelty_status` therefore remain `NOT_ASSESSED`, and no candidate receives Tier 3 or Tier 4.

## Phase 2 outputs to use

Use `results/phase2/candidates_for_perturbation_validation.tsv` as the machine-readable handoff, and `results/phase2/TOP3_candidates.md` plus `results/phase2/PHASE1_FINAL_REPORT.md` for experimental rationale. Treat `results/novel_candidates/*` as discovery layers, not validated novelty claims. The first experiments should test necessity and sufficiency independently, with receptor/ligand specificity controls and pre-specified T-cell activation and effector readouts.
