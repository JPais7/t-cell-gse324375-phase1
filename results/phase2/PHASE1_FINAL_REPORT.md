# Phase 1.5 final report

## 1. Executive summary
No candidate passes the closed Phase 2 gate. The leading statistical/biological hypothesis, **NK → Tnfsf4 → Tnfrsf4 → CD8**, is HOLD because RNA and ADT directions conflict. This is an association, not a causal conclusion.

## 2. Dataset and inferential unit
GSE324375 is analysed with the HTO-demultiplexed mouse as the inferential unit; cells measure sampling depth, not independent replication.

## 3. Activation/state discovery
The activation-like state was discovered without curated activation genes as clustering features. This provides feature independence, not external-dataset independence.

## 4. Animal-level replication
TOP1 uses 43 informative mice: 31 supporting and 12 opposing; direction consistency=0.721; median mouse effect=-0.0001864; IQR=0.001879.

## 5. Candidate discovery
Extracellular, multicellular, intrinsic-expression, inferred-TF and metabolic hypotheses remain explicitly distinct. No class balance or novelty bonus is imposed.

## 6. Biological entity validation
UniProtKB mouse annotations support TNFSF4 as a ligand and TNFRSF4 as its receptor. PKM→CD44 is invalid as conventional LR; ADAM10/17 require a documented substrate; ITGAV is not interpreted as a soluble ligand; THY1 is contact/adhesion.

## 7. Statistical robustness
The decomposable score separates replication, effect, statistics, RNA, ADT, state, interaction, temporal, LOOCV, null, biological validity, chain evidence, independence and testability. FDR is one non-dominant component.

## 8. True LOOCV
TOP1 remains top-10 in 100.0% of candidate-specific removal folds. Every fold removes one informative mouse and recalculates score, eligibility and competitive rank.

## 9. Permutation null
Within-candidate sign permutation preserves each mouse-effect magnitude while randomizing direction (1,000 permutations; seed 17). TOP1 empirical p=0.003996, empirical FDR=0.0149. This is not a biological or global-pathway null.

## 10. RNA/ADT evidence
TOP1 RNA=True; ADT=True; state=False. These modalities support association, not causality.

## 11. External perturbation evidence
TOP1: NOT_ASSESSED. No independent system-matched perturbation is currently integrated. GSE289772 is contextual pharmacological evidence for intrinsic PKM2 biology only.

## 12. Mechanistic chain
Source/signal=SUPPORTED; receptor=SUPPORTED; T-cell state=PARTIALLY_SUPPORTED; transcriptional program=SUPPORTED; effector response=PARTIALLY_SUPPORTED; myeloid response=MISSING; tumor killing=MISSING. Chain evidence coherence=0.571. Missing downstream steps are not filled with prior-paper knowledge.

## 13. Final ranking
One candidate passes all gates. Lower-scoring candidates remain HOLD or REJECT rather than being promoted to fill a list.

## 14. TOP3
1. NK → Tnfsf4 → Tnfrsf4 → CD8: class=multicellular; effect=0.2881; FDR=0.06182; mice=43; consistency=0.721; LOOCV=1.000; empirical p/FDR=0.003996/0.0149; validity=VALID_LR; external=NOT_ASSESSED; tier=2; decision=HOLD.
2. NOT AVAILABLE — did not pass all closed gates.
3. NOT AVAILABLE — did not pass all closed gates.

## 15. Why the candidate beat alternatives
It combines a biologically defensible TNFSF4–TNFRSF4 direction, RNA+ADT evidence, high animal replication, candidate-specific LOOCV and an interpretable null. Intrinsic TF hypotheses fell below the revised LOOCV gate; invalid/unknown LR interpretations were excluded.

## 16. Experimental validation
Necessity: independent TNFSF4 loss/blockade in NK cells and TNFRSF4 CRISPRi in CD8 cells. Sufficiency: physiological cross-linked TNFSF4 or TNFSF4-high NK cells. Rescue: restore TNFSF4 or use OX40 agonism, requiring TNFRSF4. Measure CD69/CD137, state program, IFNG/TNF, viability, proliferation, cytotoxicity and live tumor killing.

## 17. Falsifiers
Reject if verified independent perturbations have no concordant target-specific effect, if rescue fails, or if effects disappear after controlling viable cell counts, proliferation and generic stress.

## 18. Limitations
Cell abundance, stress, proliferation, batch, treatment and library composition remain possible residual explanations. The animal-level design and stratification reduce but cannot eliminate them. Spatial contact, longitudinal same-animal data and direct perturbation are missing.

## 19. Reproducibility
Numbered scripts generate the annotation, audit, LOOCV, null, ranking, readiness and QC outputs. UniProtKB annotations are cached; seed=17; permutations=1,000.

## 20. Final decision
Q1 biological validity: YES. Q2 mouse-removal survival: YES. Q3 mouse-level null: YES. Q4 molecular interpretation: defensible LR hypothesis. Q5 direction: NOT YET—RNA/ADT discordant. Q6 independent modalities: measured but discordant. Q7 independent perturbation: NOT_ASSESSED. Q8 testable: YES. Q9 clear falsifier: YES. Q10 residual confounding: possible. Do not begin Phase 2 until direction is resolved.
