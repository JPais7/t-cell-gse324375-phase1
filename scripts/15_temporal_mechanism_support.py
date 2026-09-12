#!/usr/bin/env python3
"""Add temporal-order support to ligand-receptor hypotheses without claiming causality."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
MECH = ROOT / "results/mechanisms/all_ranked_mechanisms.tsv.gz"
DE = ROOT / "results/experimental_contrasts/mouse_level_RNA_contrasts.tsv.gz"
OUT = ROOT / "results/mechanisms"

m = pd.read_csv(MECH, sep="\t")
d = pd.read_csv(DE, sep="\t")

# Early source expression cannot be linked longitudinally to the same mouse here.
# We therefore use a strict ordering criterion: receptor differential at 12 h and
# downstream effector/cytotoxic genes at 48 h, both relevant vs irrelevant.
# The meta-analysis combines time strata, so derive time-specific evidence from the
# explicit per-stratum DE in script 14's next-pass output when available.
specific_path = ROOT / "results/experimental_contrasts/stratum_RNA_contrasts.tsv.gz"
if specific_path.exists():
    s = pd.read_csv(specific_path, sep="\t")
else:
    s = pd.DataFrame()

if len(s):
    e12 = s[(s.contrast == "relevant_vs_irrelevant") & (s.time_hours == 12)]
    e48 = s[(s.contrast == "relevant_vs_irrelevant") & (s.time_hours == 48)]
    e12 = e12.groupby(["lineage", "feature"], observed=True).agg(early_effect=("effect_a_minus_b", "mean"), early_FDR=("FDR", "max")).reset_index()
    e48 = e48.groupby(["lineage", "feature"], observed=True).agg(late_effect=("effect_a_minus_b", "mean"), late_FDR=("FDR", "max")).reset_index()
    m = m.merge(e12, left_on=["target_lineage", "receptor"], right_on=["lineage", "feature"], how="left").drop(columns=["lineage", "feature"])
    effectors = ["Ifng", "Tnf", "Csf2", "Gzmb", "Prf1", "Fasl", "Ccl3", "Ccl4", "Xcl1"]
    late = e48[e48.feature.isin(effectors)].groupby("lineage", observed=True).agg(
        late_effector_count=("feature", lambda x: x.nunique()),
        late_effector_positive=("late_effect", lambda x: int((x > 0).sum())),
        late_effector_significant_positive=("late_effect", lambda x: int(((x > 0) & (e48.loc[x.index, "late_FDR"] < .1)).sum())),
        late_effector_best_FDR=("late_FDR", "min"),
    ).reset_index()
    m = m.merge(late, left_on="target_lineage", right_on="lineage", how="left").drop(columns="lineage")
    m["temporal_support"] = (m.early_effect > 0) & (m.early_FDR < .1) & (m.late_effector_significant_positive >= 2)
else:
    m["temporal_support"] = False
    m["temporal_note"] = "Unavailable: no same-animal longitudinal sampling; do not infer 12→48 h causality."

m["temporal_causality_label"] = "ordered cross-sectional support only"
m.to_csv(OUT / "mechanisms_with_temporal_support.tsv.gz", sep="\t", index=False)
print(m.temporal_support.value_counts(dropna=False).to_string())
