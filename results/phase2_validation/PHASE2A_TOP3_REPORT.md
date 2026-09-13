# Phase 2A TOP3 computational validation

## 1. Objective
## 2. Frozen Phase 1 state
## 3. Dataset and inferential unit
## 4. Candidate 1 — NK → Spp1 → S1pr1 → CD8
## 5. Candidate 2 — T_cell → Cd28 → Cd86 → CD8
## 6. Candidate 3 — T_cell → Lamc1 → Itga2_Itgb1 → CD4
## 7. Mouse-level associations
## 8. Adjusted models
## 9. Composition confounding
## 10. Treatment adjustment
## 11. Condition robustness
## 12. Generic activation confounding
## 13. Leave-one-mouse-out robustness
## 14. Permutation null
## 15. Multimodal evidence
## 16. Biological plausibility audit
## 17. Comparative evidence matrix
## 18. Experimental prioritization
## 19. Limitations
## 20. Phase 2A decision

C1: BACKUP_EXPERIMENTAL_PRIORITY. Raw Spearman rho -0.588, FDR 3.51e-05, bootstrap CI [-0.760,-0.343], residualized rho -0.390 (p=0.007), robust LOMO. Higher candidate-axis values are inversely associated with CD8 activation; this is not evidence of direct SPP1-S1PR1 activation.

C2: DO_NOT_PRIORITIZE_YET. Raw rho -0.367 (FDR 0.011), but treatment/condition adjustment attenuates the signal and residualized rho is -0.072 (p=0.629), consistent with between-condition structure.

C3: NOT_ASSESSABLE. 0 complete candidate-axis cases; no inferential Phase 2A analysis was performed.

PHASE2A_STATISTICAL_VALIDATION_COMPLETE = YES
