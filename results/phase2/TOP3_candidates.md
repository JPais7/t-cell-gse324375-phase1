# TOP 3 candidates

## 1. fibroblast → Adam10 → Cd44 → CD8

- Mechanism: cell-cell / extracellular signal; source=fibroblast; target=CD8; family=LR|fibroblast|Adam10.
- Evidence: ADT-supported; effect=0.321012842388682; FDR=0.0113532028462642; mice=34.0; power=HIGH; temporal=False; LOOCV=0.9705882352941176; null percentile=1.0.
- Interpretation: associated with activation, not causal.
- Necessity/sufficiency: ligand knockdown or receptor CRISPRi with rescue. Falsifier: no pre-specified activation/effector change after on-target perturbation.
- Weakness: observational and incomplete spatial/modality coverage.

## 2. Stat4

- Mechanism: inferred transcription-factor activity; source=T cell; target=CD8; family=INTRINSIC|Stat4.
- Evidence: RNA-supported; effect=0.6320756352928133; FDR=3.934714292939858e-06; mice=nan; power=MODERATE; temporal=False; LOOCV=nan; null percentile=0.397.
- Interpretation: associated with activation, not causal.
- Necessity/sufficiency: gene knockout/CRISPRi and CRISPRa rescue. Falsifier: no pre-specified activation/effector change after on-target perturbation.
- Weakness: observational and incomplete spatial/modality coverage.

## 3. T_cell → Pkm → Cd44 → CD8

- Mechanism: source-cell ↔ T-cell program / feedback; source=T_cell; target=CD8; family=MULTI|T_cell|Pkm.
- Evidence: RNA+ADT-supported; effect=0.4828625014947049; FDR=0.0394332468179697; mice=47.0; power=HIGH; temporal=False; LOOCV=0.9787234042553191; null percentile=1.0.
- Interpretation: associated with activation, not causal.
- Necessity/sufficiency: gene knockout/CRISPRi and CRISPRa rescue. Falsifier: no pre-specified activation/effector change after on-target perturbation.
- Weakness: observational and incomplete spatial/modality coverage.

