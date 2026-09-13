# Phase 1 final report

## 1. Executive summary
GSE324375 supports animal-aware, associative hypotheses for intratumoral T-cell activation. No result from this dataset alone is causal.

## 2–4. Dataset, biological replicate and QC
The HTO-demultiplexed mouse is the inferential unit. Cell counts describe sampling depth, not replication. Reversible QC, doublet calls and RNA/ADT matrices are retained.

## 5–6. Data-driven state discovery and animal validation
Unsupervised states were discovered after excluding curated activation features. The activation-like label is an external biological interpretation. State proportions and activation summaries are reported per mouse.

## 7–8. Candidate discovery and family collapsing
Candidates are separated into extracellular, intrinsic and multicellular entities. Interaction families share source and ligand; intrinsic families retain the regulator identity. No additional pathway family is invented without evidence.

## 9–11. Statistical validation, true LOOCV and null model
Priority candidates undergo actual mouse removal, score recalculation and reranking. LOOCV stability is the fraction of removal runs remaining in the top 10. The null uses within-candidate mouse-effect sign permutations (1,000; seed 17); empirical FDR remains separate from molecular FDR.

## 12–14. RNA/ADT, external perturbation and mechanistic coherence
RNA and ADT are independent components. External perturbation is counted only for verified interventions and context is retained. Each mechanistic chain records supported, partial, missing or non-applicable steps.

## 15–17. Final ranking and causality tiers
The score averages replication, effect, molecular statistics, RNA, ADT, state, temporal, interaction, LOOCV, null, coherence and curated-score independence. Verified external perturbation is added without penalizing unassessed candidates. Tier 1 is observational, Tier 2 multi-layer, Tier 3 external perturbation, and Tier 4 direct causality in the relevant system.

## 18. Phase 2 experiments
Test necessity, sufficiency and rescue using candidate-specific source/T-cell perturbations, non-targeting controls, specificity controls and pre-specified activation, effector and killing readouts.

## 19–20. Limitations and reproducibility
No spatial contact or same-animal longitudinal trajectory is available. Several contrasts are underpowered. Random seed is 17; the ranking records all components and robustness outputs.

## 21. Falsification
The leading hypothesis is falsified if verified on-target perturbation and rescue do not change the pre-specified T-cell activation/effector phenotype in the expected direction.

## TOP 10

17. T_cell → Pkm → Cd44 → CD8 (multicellular; score=0.741; power=HIGH)
910. NK → Tnfsf4 → Tnfrsf4 → CD8 (multicellular; score=0.739; power=HIGH)
770. T_cell → Adam17 → Itgb1 → CD8 (multicellular; score=0.644; power=HIGH)
2821. neutrophil → Gpi1 → Amfr → CD8 (extracellular; score=0.632; power=HIGH)
583. NK → Itgav → Thy1 → CD8 (multicellular; score=0.626; power=HIGH)
64971. Stat4 (intrinsic; score=0.601; power=MODERATE)
65002. Maff (intrinsic; score=0.569; power=MODERATE)
64979. Nfkb2 (intrinsic; score=0.549; power=MODERATE)
64995. Relb (intrinsic; score=0.548; power=MODERATE)
65030. Jund (intrinsic; score=0.545; power=MODERATE)

## Final decision fields
A. Strongest mechanistic hypothesis: T_cell → Pkm → Cd44 → CD8.
B. Strongest independent intrinsic hypothesis: Stat4.
C. Strongest independent multicellular/extracellular hypothesis: T_cell → Pkm → Cd44 → CD8.
D. Confidence: moderate; discovery-stage.
E. Evidence gaps: spatial, longitudinal and system-matched perturbation.
F. First experiment: necessity test for the leading mechanism.
G. Second experiment: sufficiency plus rescue.
H. Strongest falsifier: no effect after verified perturbation and rescue.
I. Unknown: causal direction, contact, dose and therapeutic relevance.
