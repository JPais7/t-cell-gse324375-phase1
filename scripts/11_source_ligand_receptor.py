#!/usr/bin/env python3
"""Animal-level source-ligand/receptor prioritization using mouse-consensus LR pairs."""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import liana as li
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "results" / "atlas" / "GSE324375.atlas_v1.h5ad"
TCELLS = ROOT / "results" / "tcells" / "GSE324375.Tcells.refined_v1.h5ad"
ACTIVATION = ROOT / "results" / "activation"
OUT = ROOT / "results" / "interactions"
OUT.mkdir(parents=True, exist_ok=True)
MIN_SOURCE_CELLS = 20


def components(complex_name: str) -> list[str]:
    return complex_name.split("_")


def bh(values: pd.Series) -> np.ndarray:
    arr = values.to_numpy(float)
    valid = np.isfinite(arr)
    result = np.full(len(arr), np.nan)
    result[valid] = multipletests(arr[valid], method="fdr_bh")[1]
    return result


atlas = ad.read_h5ad(ATLAS)
tcells = ad.read_h5ad(TCELLS)
resource = li.resource.select_resource("mouseconsensus").drop_duplicates()

# Retain interactions whose complete ligand and receptor complexes are measured.
atlas_genes = set(atlas.var_names)
tcell_genes = set(tcells.var_names)
resource = resource[
    resource["ligand"].map(lambda x: set(components(x)) <= atlas_genes)
    & resource["receptor"].map(lambda x: set(components(x)) <= tcell_genes)
].copy()
resource.to_csv(OUT / "mouseconsensus_measured_interactions.tsv", sep="\t", index=False)

ligand_components = sorted({g for ligand in resource.ligand for g in components(ligand)})
ligand_index = atlas.var_names.get_indexer(ligand_components)
expression = sparse.csr_matrix(atlas.X[:, ligand_index])
detection = expression.copy()
detection.data = np.ones_like(detection.data)

# Aggregate source-cell expression by mouse and broad population.
source_obs = atlas.obs[["mouse_id", "provisional_cell_type", "treatment", "time_hours", "checkpoint_blockade"]].copy()
source_obs["source_group"] = source_obs["mouse_id"].astype(str) + "|" + source_obs["provisional_cell_type"].astype(str)
groups = pd.Index(pd.unique(source_obs["source_group"]))
codes = pd.Categorical(source_obs["source_group"], categories=groups).codes
membership = sparse.csr_matrix(
    (np.ones(len(codes)), (codes, np.arange(len(codes)))), shape=(len(groups), len(codes))
)
cell_counts = np.asarray(membership.sum(axis=1)).ravel()
means = (membership @ expression).toarray() / cell_counts[:, None]
fractions = (membership @ detection).toarray() / cell_counts[:, None]

group_meta = (
    source_obs.groupby("source_group", observed=True)
    .agg(
        mouse_id=("mouse_id", "first"), source_cell_type=("provisional_cell_type", "first"),
        treatment=("treatment", "first"), time_hours=("time_hours", "first"),
        checkpoint_blockade=("checkpoint_blockade", "first"),
    )
    .loc[groups]
)
group_meta["n_source_cells"] = cell_counts
condition = (
    group_meta["treatment"].astype(str) + "|" + group_meta["time_hours"].astype(str)
    + "|" + group_meta["checkpoint_blockade"].astype(str)
)

# T-cell activation and receptor prevalence by animal and lineage.
activation_by_mouse = {}
receptor_prevalence = {}
t_counts = sparse.csr_matrix(tcells.layers["counts"])
for lineage in ("CD4", "CD8"):
    mask = tcells.obs["t_lineage_provisional"].eq(lineage).to_numpy()
    scores = tcells.obs.loc[mask].groupby("mouse_id", observed=True)[
        ["score_activation_immediate", "score_activation_ADT"]
    ].median()
    activation_by_mouse[lineage] = scores
    lineage_counts = t_counts[mask]
    receptor_prevalence[lineage] = {
        gene: float((lineage_counts[:, tcells.var_names.get_loc(gene)] > 0).mean())
        for gene in sorted({g for receptor in resource.receptor for g in components(receptor)})
    }

# Mouse-level correlation of source ligand expression with T-cell activation.
association_rows = []
gene_position = {gene: i for i, gene in enumerate(ligand_components)}
for source_type in sorted(group_meta.source_cell_type.unique()):
    source_mask = (group_meta.source_cell_type == source_type) & (group_meta.n_source_cells >= MIN_SOURCE_CELLS)
    source_indices = np.flatnonzero(source_mask.to_numpy())
    if len(source_indices) < 6:
        continue
    source_mice = group_meta.iloc[source_indices].mouse_id.astype(str)
    source_condition = condition.iloc[source_indices]
    for ligand in sorted(resource.ligand.unique()):
        positions = [gene_position[g] for g in components(ligand)]
        ligand_mean = means[source_indices][:, positions].min(axis=1)
        ligand_fraction = fractions[source_indices][:, positions].min(axis=1)
        for lineage in ("CD4", "CD8"):
            scores = activation_by_mouse[lineage]
            common = [i for i, mouse in enumerate(source_mice) if mouse in scores.index]
            if len(common) < 6:
                continue
            mice = source_mice.iloc[common].tolist()
            x = ligand_mean[common]
            y_rna = scores.loc[mice, "score_activation_immediate"].to_numpy()
            y_adt = scores.loc[mice, "score_activation_ADT"].to_numpy()
            strata = source_condition.iloc[common].to_numpy()
            # Residualize within experimental strata so treatment is not the only
            # source of an apparent ligand–activation relationship.
            x_res = x - pd.Series(x).groupby(strata).transform("mean").to_numpy()
            rna_res = y_rna - pd.Series(y_rna).groupby(strata).transform("mean").to_numpy()
            adt_res = y_adt - pd.Series(y_adt).groupby(strata).transform("mean").to_numpy()
            rho_rna, p_rna = (spearmanr(x_res, rna_res) if np.std(x_res) > 0 and np.std(rna_res) > 0 else (np.nan, np.nan))
            rho_adt, p_adt = (spearmanr(x_res, adt_res) if np.std(x_res) > 0 and np.std(adt_res) > 0 else (np.nan, np.nan))
            association_rows.append(
                {
                    "source_cell_type": source_type, "ligand": ligand, "target_lineage": lineage,
                    "n_mice": len(common), "source_mean_expression": float(np.mean(x)),
                    "source_detection_fraction": float(np.mean(ligand_fraction[common])),
                    "rho_condition_adjusted_RNA_activation": rho_rna, "p_RNA": p_rna,
                    "rho_condition_adjusted_ADT_activation": rho_adt, "p_ADT": p_adt,
                }
            )

associations = pd.DataFrame(association_rows)
associations["FDR_RNA"] = associations.groupby("target_lineage", observed=True)["p_RNA"].transform(bh)
associations["FDR_ADT"] = associations.groupby("target_lineage", observed=True)["p_ADT"].transform(bh)
associations.to_csv(OUT / "source_ligand_mouse_associations.tsv.gz", sep="\t", index=False)

# Join source evidence to receptor activation effects.
rankings = pd.read_csv(ACTIVATION / "cross_modal_activation_gene_ranking.tsv.gz", sep="\t")
interaction_rows = []
for row in resource.itertuples(index=False):
    ligand_assoc = associations[associations.ligand == row.ligand]
    for lineage in ("CD4", "CD8"):
        receptor_genes = components(row.receptor)
        receptor_rank = rankings[(rankings.lineage == lineage) & rankings.gene.isin(receptor_genes)]
        if len(receptor_rank) != len(receptor_genes):
            continue
        receptor_effect_rna = float(receptor_rank.median_paired_log2FC_RNA.min())
        receptor_effect_adt = float(receptor_rank.median_paired_log2FC_ADT.min())
        receptor_fdr_rna = float(receptor_rank.FDR_RNA.max())
        receptor_fdr_adt = float(receptor_rank.FDR_ADT.max())
        receptor_detect = min(receptor_prevalence[lineage][g] for g in receptor_genes)
        for assoc in ligand_assoc[ligand_assoc.target_lineage == lineage].itertuples(index=False):
            cross_modal_direction = (
                assoc.rho_condition_adjusted_RNA_activation > 0
                and assoc.rho_condition_adjusted_ADT_activation > 0
                and receptor_effect_rna > 0 and receptor_effect_adt > 0
            )
            interaction_rows.append(
                {
                    "source_cell_type": assoc.source_cell_type, "ligand": row.ligand,
                    "receptor": row.receptor, "target_lineage": lineage, "n_mice": assoc.n_mice,
                    "source_detection_fraction": assoc.source_detection_fraction,
                    "target_receptor_detection_fraction": receptor_detect,
                    "ligand_rho_RNA": assoc.rho_condition_adjusted_RNA_activation,
                    "ligand_FDR_RNA": assoc.FDR_RNA,
                    "ligand_rho_ADT": assoc.rho_condition_adjusted_ADT_activation,
                    "ligand_FDR_ADT": assoc.FDR_ADT,
                    "receptor_log2FC_RNA_selector": receptor_effect_rna,
                    "receptor_FDR_RNA_selector": receptor_fdr_rna,
                    "receptor_log2FC_ADT_selector": receptor_effect_adt,
                    "receptor_FDR_ADT_selector": receptor_fdr_adt,
                    "cross_modal_positive_direction": cross_modal_direction,
                    "evidence_class": "C — correlational",
                }
            )

interactions = pd.DataFrame(interaction_rows)
interactions["association_component"] = interactions[["ligand_rho_RNA", "ligand_rho_ADT"]].clip(lower=0).mean(axis=1)
interactions["receptor_component"] = interactions[
    ["receptor_log2FC_RNA_selector", "receptor_log2FC_ADT_selector"]
].clip(lower=0).mean(axis=1).clip(upper=3) / 3
interactions["expression_component"] = np.sqrt(
    interactions.source_detection_fraction * interactions.target_receptor_detection_fraction
)
interactions["replication_component"] = (
    (interactions.ligand_FDR_RNA < 0.1).astype(float)
    + (interactions.ligand_FDR_ADT < 0.1).astype(float)
    + (interactions.receptor_FDR_RNA_selector < 0.05).astype(float)
    + (interactions.receptor_FDR_ADT_selector < 0.05).astype(float)
) / 4
interactions["interaction_score"] = (
    0.35 * interactions.association_component
    + 0.25 * interactions.receptor_component
    + 0.20 * interactions.expression_component
    + 0.20 * interactions.replication_component
)
basic_support = (
    interactions.cross_modal_positive_direction
    & (interactions.receptor_FDR_RNA_selector < 0.05)
    & (interactions.receptor_FDR_ADT_selector < 0.05)
    & (interactions.source_detection_fraction >= 0.05)
    & (interactions.target_receptor_detection_fraction >= 0.05)
)
both_ligand_modalities = (interactions.ligand_FDR_RNA < 0.1) & (interactions.ligand_FDR_ADT < 0.1)
one_ligand_modality = (interactions.ligand_FDR_RNA < 0.1) | (interactions.ligand_FDR_ADT < 0.1)
interactions["evidence_tier"] = np.select(
    [basic_support & both_ligand_modalities, basic_support & one_ligand_modality],
    ["Tier 1 — replicated multimodal", "Tier 2 — single-modality source support"],
    default="Tier 3 — exploratory",
)
interactions["tier_rank"] = interactions.evidence_tier.map(
    {"Tier 1 — replicated multimodal": 1, "Tier 2 — single-modality source support": 2,
     "Tier 3 — exploratory": 3}
)
interactions.sort_values(["tier_rank", "interaction_score"], ascending=[True, False], inplace=True)
interactions.to_csv(OUT / "ranked_source_ligand_receptor_candidates.tsv.gz", sep="\t", index=False)
interactions.head(200).to_csv(OUT / "top200_source_ligand_receptor_candidates.tsv", sep="\t", index=False)
interactions[interactions.tier_rank <= 2].to_csv(
    OUT / "supported_source_ligand_receptor_candidates.tsv", sep="\t", index=False
)
print(f"Measured mouse-consensus interactions: {len(resource):,}")
print(f"Ranked source→T-cell hypotheses: {len(interactions):,}")
