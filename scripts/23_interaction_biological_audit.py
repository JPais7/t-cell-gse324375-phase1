#!/usr/bin/env python3
"""Apply UniProt-derived entity annotations to every proposed interaction."""
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/phase2"
x=pd.read_csv(OUT/"candidate_evidence_pre_audit.tsv",sep="\t",low_memory=False); lr=x[x.candidate_type.isin(["extracellular","multicellular"])].copy()
ann=pd.read_csv(OUT/"interaction_entity_annotation.tsv",sep="\t").set_index("gene")
def entities(value): return [p for p in str(value).split("_") if p and p!="nan"]
def rec(value): return ann.reindex(entities(value))
def join(value,col): return " | ".join(rec(value)[col].fillna("unknown/unvalidated" if col=="entity_type" else "").astype(str))
def allflag(value,col):
 z=rec(value); return bool(len(z) and z[col].fillna(False).astype(bool).all())
def anyflag(value,col): return bool(rec(value)[col].fillna(False).astype(bool).any())
lr["ligand_entity_type"]=lr.ligand.map(lambda g:join(g,"entity_type")); lr["receptor_entity_type"]=lr.receptor.map(lambda g:join(g,"entity_type")); lr["ligand_localization"]=lr.ligand.map(lambda g:join(g,"cellular_localization")); lr["receptor_localization"]=lr.receptor.map(lambda g:join(g,"cellular_localization"))
lr["ligand_supported_annotation"]=lr.ligand.map(lambda g:allflag(g,"ligand_supported")); lr["receptor_supported_annotation"]=lr.receptor.map(lambda g:allflag(g,"receptor_supported"))
intrinsic_lig=lr.ligand.map(lambda g:anyflag(g,"metabolic_regulator") or anyflag(g,"transcription_factor") or (allflag(g,"enzyme") and not anyflag(g,"secreted")))
inverted=lr.ligand.map(lambda g:all(t=="receptor" for t in rec(g).entity_type.fillna("unknown/unvalidated")))
shedding=lr.ligand.map(lambda g:anyflag(g,"shedding_or_processing_factor")); contact=lr.ligand.map(lambda g:anyflag(g,"adhesion_molecule") and (anyflag(g,"membrane_associated") or anyflag(g,"secreted"))); canonical=lr.ligand_supported_annotation&lr.receptor_supported_annotation
def substrate_supported(row):
 notes=" ".join(rec(row.ligand).biological_notes.fillna("").astype(str)).lower(); return any(p.lower() in notes for p in entities(row.receptor))
substrate=lr.apply(substrate_supported,axis=1)
lr["interaction_validity"]=np.select([intrinsic_lig|inverted,shedding&substrate,shedding,canonical,contact&lr.receptor_supported_annotation],["INVALID_LR","SHEDDING_PROCESSING","UNKNOWN","VALID_LR","MEMBRANE_CONTACT"],default="UNKNOWN")
lr["interaction_type"]=lr.interaction_validity; lr["direction_validity"]=np.select([intrinsic_lig|inverted,shedding&substrate,shedding,canonical,contact&lr.receptor_supported_annotation],["INVALID","PLAUSIBLE_CONTACT_OR_PROCESSING","UNKNOWN","SUPPORTED","PLAUSIBLE_CONTACT_OR_PROCESSING"],default="UNKNOWN")
lr["validation_source"]="UniProtKB mouse entity annotation plus interaction-resource hypothesis"; lr["validation_confidence"]=np.select([lr.interaction_validity.isin(["INVALID_LR","VALID_LR"]),lr.interaction_validity.isin(["MEMBRANE_CONTACT","SHEDDING_PROCESSING"])],["HIGH","MODERATE"],default="LOW")
lr["biological_reinterpretation"]=np.select([lr.ligand.eq("Pkm"),lr.ligand.isin(["Adam10","Adam17"]),lr.ligand.eq("Itgav"),lr.ligand.eq("Thy1")],["T-cell intrinsic PKM-associated metabolic state; extracellular PKM→CD44 is not established","Membrane protease/shedding hypothesis; retain only with an annotated substrate","ITGAV is an integrin receptor/contact molecule and is directionally inappropriate as a conventional soluble ligand","THY1 is a membrane adhesion/contact molecule; interpret only as contact-dependent"],default="")
lr["valid_for_mechanistic_ranking"]=lr.interaction_validity.isin(["VALID_LR","MEMBRANE_CONTACT","SHEDDING_PROCESSING"])
audit_cols=["candidate","candidate_type","source","target","ligand","receptor","ligand_entity_type","receptor_entity_type","ligand_localization","receptor_localization","ligand_supported_annotation","receptor_supported_annotation","interaction_type","interaction_validity","direction_validity","validation_source","validation_confidence","biological_reinterpretation","valid_for_mechanistic_ranking"]
lr[audit_cols].to_csv(OUT/"interaction_biological_audit.tsv",sep="\t",index=False)
(OUT/"interaction_biological_audit.md").write_text("# Interaction biological audit\n\nEntity types and localizations are derived reproducibly from the cached UniProtKB Mus musculus export. Interaction-resource membership is treated only as a hypothesis. Valid conventional LR, membrane contact, and substrate-supported shedding mechanisms are distinguished from invalid and unknown pairs. PKM→CD44 is excluded as LR; ADAM10/17 require an annotated substrate; ITGAV in ligand position is not treated as soluble signaling; THY1 is contact/adhesion.\n")
print(lr.interaction_validity.value_counts().to_string())
