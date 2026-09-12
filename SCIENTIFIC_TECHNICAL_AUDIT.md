# Scientific and technical audit — GSE324375 Phase 1

## What is sound

- HTO-demultiplexed mouse/library-hash units are retained as the biological replicate.
- RNA and ADT are kept as separate activation axes; defining features are excluded from the corresponding discovery contrasts.
- QC is reversible, raw counts are preserved, and HTO and transcriptomic doublet evidence are recorded independently.
- Source→ligand→receptor analyses are condition-adjusted, expression-aware and explicitly correlational.
- Experimental contrasts are evaluated within comparable time/checkpoint strata rather than pooling incompatible conditions.

## Issues corrected in this review

1. Added an unsupervised T-cell atlas built after removing curated activation, effector, proliferation, dysfunction and interferon genes. Curated scores are used only afterward to interpret which unsupervised state is activation-enriched.
2. Added an explicit gene-level circularity table (`used_in_curated_score`, data-driven association and independent-candidate flag).
3. Added separate extracellular, intrinsic and multicellular candidate rankings, with transparent evidence columns and a novelty indicator that means non-trivial hypothesis, not absence from the dataset.
4. Added mouse-level contrast, TF-activity and sensitivity outputs. Every contrast records included strata, animal counts and a supported/fragile design label.

## Remaining limitations

- GSE324375 is cross-sectional across time; 12→48 h ordering is not longitudinal causality.
- Ligand–receptor compatibility does not establish spatial contact or signaling.
- DoRothEA/ULM TF activities are inferred regulatory scores, not direct TF measurements.
- Some relevant-vs-innate strata contain too few animals and remain LOW_POWER/fragile.
- The final candidate lists are discovery hypotheses. Perturbation is required for necessity and sufficiency.
