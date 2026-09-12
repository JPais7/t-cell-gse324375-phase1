# GSE324375 — Phase 1 discovery analysis

Goal: identify activation-associated signals in intratumoral T cells while
preserving the distinction between association, experimental evidence and
causality.

## Analysis principles

- The biological replicate is the HTO-demultiplexed mouse, not a cell or GEO file.
- QC thresholds will be reported and justified; unfiltered and filtered objects
  will both be retained.
- Differential expression will use mouse-level pseudobulk where replication allows.
- Activation scores will be validated across animals/conditions, and score genes
  will not be rediscovered through circular differential-expression tests.
- Ligand–receptor results are interaction hypotheses unless supported by an
  experimental contrast and downstream response.

## Layout

- `data/raw/`: downloaded GEO archive (not for version control)
- `metadata/`: GEO MINiML, file manifest, antibody/HTO feature definitions
- `scripts/`: reproducible numbered analysis scripts
- `results/qc/`: Phase 1 QC tables and figures
