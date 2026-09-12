#!/usr/bin/env python3
"""Combine animal-level and LIANA evidence into an auditable mechanism shortlist."""

from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
IN = ROOT / "results" / "interactions"
OUT = ROOT / "results" / "mechanisms"
OUT.mkdir(parents=True, exist_ok=True)

animal = pd.read_csv(IN / "ranked_source_ligand_receptor_candidates.tsv.gz", sep="\t")
liana = pd.read_csv(IN / "liana_sources_to_Tcells.tsv.gz", sep="\t").rename(
    columns={
        "source": "source_cell_type",
        "ligand_complex": "ligand",
        "receptor_complex": "receptor",
        "magnitude_rank": "liana_magnitude_rank",
        "specificity_rank": "liana_specificity_rank",
    }
)
liana = liana.sort_values("liana_magnitude_rank").drop_duplicates(
    ["source_cell_type", "ligand", "receptor"]
)
merged = animal.merge(
    liana[
        ["source_cell_type", "ligand", "receptor", "liana_magnitude_rank", "liana_specificity_rank"]
    ],
    on=["source_cell_type", "ligand", "receptor"],
    how="left",
)
merged["liana_component"] = (1 - merged["liana_magnitude_rank"].fillna(1)).clip(0, 1)
merged["final_discovery_score"] = 0.80 * merged["interaction_score"] + 0.20 * merged["liana_component"]
merged["mechanistic_chain"] = (
    merged.source_cell_type.astype(str) + " → " + merged.ligand.astype(str) + " → "
    + merged.receptor.astype(str) + " → " + merged.target_lineage.astype(str)
    + " activation-associated program"
)
merged["causality"] = "CORRELATIONAL"
merged["necessity_test"] = (
    "Perturb ligand in " + merged.source_cell_type.astype(str)
    + " and independently KO/CRISPRi receptor in target T cells"
)
merged["sufficiency_test"] = (
    "Increase ligand exposure or expression, with receptor-deficient T cells as specificity control"
)
merged["score_definition"] = (
    "0.80×animal-level interaction score + 0.20×(1−LIANA magnitude rank)"
)

# Preserve every candidate, then generate a deliberately mixed shortlist: all
# supported Tier-2 interactions followed by expression-plausible exploratory ones.
merged.sort_values(["tier_rank", "final_discovery_score"], ascending=[True, False], inplace=True)
merged.to_csv(OUT / "all_ranked_mechanisms.tsv.gz", sep="\t", index=False)

supported = merged[merged.tier_rank <= 2]
exploratory = merged[
    (merged.tier_rank == 3)
    & (merged.liana_magnitude_rank < 0.05)
    & (merged.source_detection_fraction >= 0.10)
    & (merged.target_receptor_detection_fraction >= 0.10)
    & ((merged.ligand_rho_RNA > 0) | (merged.ligand_rho_ADT > 0))
    & (merged.receptor_log2FC_RNA_selector > 0)
    & (merged.receptor_log2FC_ADT_selector > 0)
].sort_values("final_discovery_score", ascending=False)

shortlist = pd.concat([supported, exploratory], ignore_index=True).drop_duplicates(
    ["source_cell_type", "ligand", "receptor", "target_lineage"]
).head(20)
shortlist.insert(0, "candidate_rank", np.arange(1, len(shortlist) + 1))
shortlist.to_csv(OUT / "mechanism_shortlist_20.tsv", sep="\t", index=False)

summary_columns = [
    "candidate_rank", "evidence_tier", "mechanistic_chain", "final_discovery_score",
    "ligand_rho_RNA", "ligand_FDR_RNA", "ligand_rho_ADT", "ligand_FDR_ADT",
    "receptor_log2FC_RNA_selector", "receptor_FDR_RNA_selector",
    "receptor_log2FC_ADT_selector", "receptor_FDR_ADT_selector",
    "liana_magnitude_rank", "causality", "necessity_test", "sufficiency_test",
]
shortlist[summary_columns].to_csv(OUT / "mechanism_shortlist_20_compact.tsv", sep="\t", index=False)
print(shortlist[summary_columns[:5]].to_string(index=False))
