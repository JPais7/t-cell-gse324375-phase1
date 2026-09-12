#!/usr/bin/env python3
"""Expression-based LIANA consensus; complementary, not causal evidence."""

from pathlib import Path

import anndata as ad
import liana as li


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results" / "atlas" / "GSE324375.atlas_v1.h5ad"
OUT = ROOT / "results" / "interactions"

adata = ad.read_h5ad(SOURCE)
result = li.mt.rank_aggregate(
    adata,
    groupby="provisional_cell_type",
    resource_name="mouseconsensus",
    expr_prop=0.10,
    min_cells=20,
    use_raw=False,
    n_perms=1000,
    n_jobs=4,
    inplace=False,
    verbose=True,
)
result["evidence_class"] = "C — expression plausibility"
result.to_csv(OUT / "liana_global_expression_plausibility.tsv.gz", sep="\t", index=False)
target = result[result["target"] == "T_cell"].copy()
target.to_csv(OUT / "liana_sources_to_Tcells.tsv.gz", sep="\t", index=False)
print(f"LIANA interactions: {len(result):,}; source→T cell: {len(target):,}")
