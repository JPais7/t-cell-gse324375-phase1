# Phase 1.5 final report

## 1. Executive summary
No candidate passes the closed Phase 2 gate. The leading intrinsic hypothesis, **Eif2ak3**, is HOLD because one or more biological, direction, robustness or causal-validation gates remain unresolved. This is an association, not a causal conclusion.

## 2. Dataset and inferential unit
GSE324375 is analysed with the HTO-demultiplexed mouse as the inferential unit; cells measure sampling depth, not independent replication.

## 3. Activation/state discovery
The activation-like state was discovered without curated activation genes as clustering features. This provides feature independence, not external-dataset independence.

## 4. Animal-level replication
TOP1 uses 42 informative mice: 0 supporting and 0 opposing; direction consistency=nan; median mouse effect=-0.5267; IQR=0.9302.

## 5. Candidate discovery
Extracellular, multicellular, intrinsic-expression, inferred-TF and metabolic hypotheses remain explicitly distinct. No class balance or novelty bonus is imposed.

## 6. Biological entity validation
UniProtKB mouse annotations support TNFSF4 as a ligand and TNFRSF4 as its receptor. PKM→CD44 is invalid as conventional LR; ADAM10/17 require a documented substrate; ITGAV is not interpreted as a soluble ligand; THY1 is contact/adhesion.

## 7. Statistical robustness
The decomposable score separates replication, effect, statistics, RNA, ADT, state, interaction, temporal, LOOCV, null, biological validity, chain evidence, independence and testability. FDR is one non-dominant component.

## 8. True LOOCV
TOP1 remains top-10 in 100.0% of candidate-specific removal folds. Every fold removes one informative mouse and recalculates score, eligibility and competitive rank.

## 9. Permutation null
Within-candidate sign permutation preserves each mouse-effect magnitude while randomizing direction (1,000 permutations; seed 17). TOP1 empirical p=0.000999, empirical FDR=0.005162. This is not a biological or global-pathway null.

## 10. RNA/ADT evidence
TOP1 RNA=True; ADT=True; state=True. These modalities support association, not causality.

## 11. External perturbation evidence
TOP1: NOT_ASSESSED. No independent system-matched perturbation is currently integrated. GSE289772 is contextual pharmacological evidence for intrinsic PKM2 biology only.

## 12. Mechanistic chain
Source/signal=NOT_APPLICABLE; receptor=NOT_APPLICABLE; T-cell state=SUPPORTED; transcriptional program=SUPPORTED; effector response=PARTIALLY_SUPPORTED; myeloid response=MISSING; tumor killing=MISSING. Chain evidence coherence=0.500. Missing downstream steps are not filled with prior-paper knowledge.

## 13. Final ranking
No candidate passes all gates. Lower-scoring candidates remain HOLD or REJECT rather than being promoted to fill a list.

## 14. TOP3
1. Eif2ak3: class=intrinsic; effect=nan; FDR=9.334e-30; mice=42; consistency=nan; LOOCV=1.000; empirical p/FDR=0.000999/0.005162; validity=nan; external=NOT_ASSESSED; tier=2; decision=HOLD.
2. NOT AVAILABLE — did not pass all closed gates.
3. NOT AVAILABLE — did not pass all closed gates.

## 15. Why the candidate beat alternatives
The leading candidate is ranked by its candidate-type-specific evidence, animal-level robustness, null calibration and biological validity. It remains observational and should not be interpreted as causal or therapeutically beneficial without perturbation.

## 16. Experimental validation
Necessity and sufficiency tests must match the candidate class. For intrinsic candidates, use two independent CRISPRi/KO reagents, physiological rescue, activation and functional readouts, plus viability, apoptosis, proliferation and stress controls.

## 17. Falsifiers
Reject if verified independent perturbations have no concordant target-specific effect, if rescue fails, or if effects disappear after controlling viable cell counts, proliferation and generic stress.

## 18. Limitations
Cell abundance, stress, proliferation, batch, treatment and library composition remain possible residual explanations. The animal-level design and stratification reduce but cannot eliminate them. Spatial contact, longitudinal same-animal data and direct perturbation are missing.

## 19. Reproducibility
Numbered scripts generate the annotation, audit, LOOCV, null, ranking, readiness and QC outputs. UniProtKB annotations are cached; seed=17; permutations=1,000.

## 20. Final decision
Candidate-specific gates remain unresolved; independent perturbational validation is required.
