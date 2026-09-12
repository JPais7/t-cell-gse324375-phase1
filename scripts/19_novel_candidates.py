#!/usr/bin/env python3
"""Balanced candidate rankings across extracellular, intrinsic and multicellular mechanisms."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/novel_candidates"; OUT.mkdir(parents=True,exist_ok=True)
audit=pd.read_csv(ROOT/"results/data_driven/gene_evidence_circularity_audit.tsv.gz",sep="\t")
tf=pd.read_csv(ROOT/"results/regulators/TF_activity_experimental_contrasts.tsv",sep="\t")
tf=tf.sort_values("FDR").groupby("TF",observed=True).head(1)
genes=pd.read_csv(ROOT/"results/experimental_contrasts/mouse_level_RNA_contrasts.tsv.gz",sep="\t")
genes=genes.sort_values("FDR").groupby(["feature","lineage"],observed=True).head(1)
prot=pd.read_csv(ROOT/"results/activation/paired_ADT_activation_results.tsv",sep="\t")
prot=prot.sort_values("FDR").groupby("protein",observed=True).head(1)
state=audit[audit.associated_with_data_driven_state].copy()
state["candidate_type"]="intrinsic"
state["mechanism_class"]="T-cell-intrinsic regulator/program"
state["effect_size"]=state.data_driven_log2FC
state["FDR"]=state.data_driven_FDR
state["RNA_evidence"]=state.RNA_FDR<0.05; state["ADT_evidence"]=state.ADT_FDR<0.05
state["data_driven_state_evidence"]=True; state["cell_source_evidence"]=False; state["receptor_evidence"]=False; state["pathway_evidence"]=False
state["treatment_consistency"]=np.nan; state["novelty_indicator"]=(~state.used_in_curated_score)
state["evidence_class"]=np.select([state.RNA_evidence&state.ADT_evidence,state.RNA_evidence,state.ADT_evidence],["RNA+ADT-supported","RNA-supported","ADT-supported"],default="state-supported")
state=state.rename(columns={"gene":"candidate"})
state["score"]=(-np.log10(state.FDR.clip(lower=1e-300))+np.maximum(state.effect_size,0)+state.novelty_indicator.astype(float))
intrinsic=state[["candidate","candidate_type","mechanism_class","effect_size","FDR","RNA_evidence","ADT_evidence","data_driven_state_evidence","cell_source_evidence","receptor_evidence","pathway_evidence","treatment_consistency","novelty_indicator","evidence_class","score"]]

# TFs are an independent intrinsic layer and are explicitly labeled inferred.
t=tf[(tf.FDR<0.1)].copy(); t["candidate"]=t.TF; t["candidate_type"]="intrinsic"; t["mechanism_class"]="inferred transcription-factor activity"; t["effect_size"]=t.activity_effect_a_minus_b; t["FDR"]=t.FDR; t["RNA_evidence"]=True; t["ADT_evidence"]=False; t["data_driven_state_evidence"]=False; t["cell_source_evidence"]=False; t["receptor_evidence"]=False; t["pathway_evidence"]=False; t["treatment_consistency"]=t.direction_consistency; t["novelty_indicator"]=~t.candidate.isin(set(state.candidate)); t["evidence_class"]="RNA-supported (inferred activity)"; t["score"]=-np.log10(t.FDR.clip(lower=1e-300))+np.maximum(t.effect_size,0)+t.novelty_indicator.astype(float)
intrinsic=pd.concat([intrinsic,t[intrinsic.columns]],ignore_index=True).drop_duplicates("candidate").sort_values("score",ascending=False)

# Extracellular: preserve the previous interaction evidence, but expose every component.
lr=pd.read_csv(ROOT/"results/mechanisms/mechanisms_with_temporal_support.tsv.gz",sep="\t")
lr["candidate"]=lr.source_cell_type+" → "+lr.ligand+" → "+lr.receptor+" → "+lr.target_lineage
lr["candidate_type"]="extracellular"; lr["mechanism_class"]="cell-cell / extracellular signal"; lr["effect_size"]=lr.final_discovery_score; lr["FDR"]=np.minimum(lr.ligand_FDR_RNA.fillna(1),lr.ligand_FDR_ADT.fillna(1)); lr["RNA_evidence"]=lr.ligand_FDR_RNA<0.1; lr["ADT_evidence"]=lr.ligand_FDR_ADT<0.1; lr["data_driven_state_evidence"]=False; lr["cell_source_evidence"]=lr.source_detection_fraction>=.1; lr["receptor_evidence"]=lr.target_receptor_detection_fraction>=.1; lr["pathway_evidence"]=False; lr["treatment_consistency"]=lr.cross_modal_positive_direction; lr["novelty_indicator"]=~lr.ligand.isin(set(["Ifng","Tnf","Il2","Il7","Ccl5"])); lr["evidence_class"]=np.select([lr.RNA_evidence&lr.ADT_evidence,lr.RNA_evidence,lr.ADT_evidence],["RNA+ADT-supported","RNA-supported","ADT-supported"],default="interaction-supported"); lr["score"]=lr.final_discovery_score+lr.novelty_indicator.astype(float)
extracellular=lr[["candidate","candidate_type","mechanism_class","effect_size","FDR","RNA_evidence","ADT_evidence","data_driven_state_evidence","cell_source_evidence","receptor_evidence","pathway_evidence","treatment_consistency","novelty_indicator","evidence_class","score"]].sort_values("score",ascending=False)

# Multicellular candidates are source/program combinations, not just LR edges.
multi=lr[lr.source_cell_type.isin(["T_cell","monocyte_macrophage","dendritic","NK","melanoma"])].copy(); multi["candidate_type"]="multicellular"; multi["mechanism_class"]="source-cell ↔ T-cell program / feedback"; multi["score"]=(multi.final_discovery_score+multi.treatment_consistency.fillna(False).astype(float)); multicellular=multi[extracellular.columns].sort_values("score",ascending=False)
for name,frame in [("extracellular",extracellular),("intrinsic",intrinsic),("multicellular",multicellular)]: frame.head(200).to_csv(OUT/f"{name}_candidates.tsv",sep="\t",index=False)
allc=pd.concat([extracellular,intrinsic,multicellular],ignore_index=True); allc.to_csv(OUT/"all_candidates.tsv.gz",sep="\t",index=False)
print(allc.groupby("candidate_type",observed=True).size().to_string())
