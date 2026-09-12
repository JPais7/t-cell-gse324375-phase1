#!/usr/bin/env python3
"""Rank-based pathway enrichment for animal-level activation effects."""

from pathlib import Path

import gseapy as gp
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "activation"
OUT = SOURCE / "pathways"
OUT.mkdir(parents=True, exist_ok=True)
LIBRARIES = ["Reactome_2022", "GO_Biological_Process_2023"]

all_results = []
for lineage in ("CD4", "CD8"):
    for selector in ("RNA_activation", "ADT_activation"):
        table = pd.read_csv(SOURCE / f"{lineage}.{selector}.paired_DE.tsv.gz", sep="\t")
        table = table[table["eligible_discovery_result"]].dropna(subset=["median_paired_log2FC"])
        table = table[table["median_paired_log2FC"].abs() > 1e-8].copy()
        table["rank_metric"] = table["median_paired_log2FC"] * (
            1 - np.log10(table["pvalue"].clip(lower=1e-300))
        )
        ranking = table[["gene", "rank_metric"]].drop_duplicates("gene")
        ranking = ranking.sort_values("rank_metric", ascending=False)
        for library in LIBRARIES:
            result = gp.prerank(
                rnk=ranking,
                gene_sets=library,
                organism="Mouse",
                min_size=10,
                max_size=500,
                permutation_num=1000,
                seed=17,
                threads=4,
                verbose=False,
                outdir=None,
            ).res2d
            result.insert(0, "lineage", lineage)
            result.insert(1, "selector", selector)
            result.insert(2, "gene_set_library", library)
            all_results.append(result)
            print(lineage, selector, library, len(result), flush=True)

pd.concat(all_results, ignore_index=True).to_csv(
    OUT / "activation_preranked_GSEA.tsv.gz", sep="\t", index=False
)
