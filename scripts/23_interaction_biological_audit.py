#!/usr/bin/env python3
"""Auditable biological-type screen for ligand–receptor candidate interpretation."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/"results/phase2"
x=pd.read_csv(OUT/"candidate_evidence_matrix.tsv",sep="\t"); lr=x[x.candidate_type.isin(["extracellular","multicellular"])].copy()
METABOLIC={"Pkm","Gpi1","Hdc","Hk2","Ldha","Eno1","Gapdh","Aldoa"}; INTRACELLULAR={"Arf6","Vim","Traf2","Smap1","Stat4","Nfkb2","Relb","Jund","Rbpj"}; MEMBRANE={"Adam10","Adam17","Itgav","Thy1","Cd72","Cd38","H2-Ab1"}
CANONICAL={"Tnfsf4","Pdcd1lg2","Il15","Ccl11","Ccl25","Cd80","Cd86","Spp1","Grn","Vegfa","Lgals3bp"}
def ligand_class(g):
 if g in METABOLIC:return "METABOLIC_ENZYME"
 if g in INTRACELLULAR:return "INTRACELLULAR_NOT_LIGAND"
 if g in MEMBRANE:return "MEMBRANE_ASSOCIATED_SIGNAL"
 if g in CANONICAL or str(g).startswith(("Ccl","Cxcl","Il","Tnfsf")):return "VALID_CANONICAL_LIGAND_RECEPTOR"
 if str(g).startswith(("Col","Lam","Hspg")):return "VALID_NONCANONICAL_EXTRACELLULAR"
 return "UNKNOWN_UNVALIDATED"
def receptor_class(g):
 parts=str(g).split("_")
 if any(p in INTRACELLULAR or p in METABOLIC for p in parts):return "INTRACELLULAR_NOT_RECEPTOR"
 if all(p.startswith(("Tnfrsf","Ccr","Cxcr","Il","Itga","Itgb")) or p in {"Cd44","Cd5","Cd28","Pdcd1","Lag3","Ptprs","Dpp4","Amfr","Pecam1","Thy1","Notch1","Rgmb"} for p in parts):return "VALID_CANONICAL_LIGAND_RECEPTOR"
 return "UNKNOWN_UNVALIDATED"
lr["ligand_biological_class"]=lr.ligand.map(ligand_class); lr["receptor_biological_class"]=lr.receptor.map(receptor_class)
invalid=lr.ligand_biological_class.isin(["METABOLIC_ENZYME","INTRACELLULAR_NOT_LIGAND"])|lr.receptor_biological_class.eq("INTRACELLULAR_NOT_RECEPTOR")|lr.ligand.eq("Itgav")
unvalidated=lr.ligand_biological_class.eq("UNKNOWN_UNVALIDATED")|lr.receptor_biological_class.eq("UNKNOWN_UNVALIDATED")|(lr.ligand.isin(["Adam10","Adam17"])&~lr.receptor.isin(["Notch1","Tnfrsf1b","Il6r","Sell"]))
lr["interaction_validity"]=np.where(invalid,"INVALID_AS_LIGAND_RECEPTOR",np.where(unvalidated,"UNKNOWN_UNVALIDATED","VALID_FOR_HYPOTHESIS"))
lr["interaction_type"]=np.where(lr.ligand_biological_class.eq("MEMBRANE_ASSOCIATED_SIGNAL"),"membrane-associated/contact or processing","ligand-receptor hypothesis"); lr["validation_source"]="UniProt/GO-informed explicit rule set; mouseconsensus remains interaction-resource evidence"; lr["validation_confidence"]=np.where(lr.interaction_validity.eq("INVALID_AS_LIGAND_RECEPTOR"),"HIGH",np.where(lr.interaction_validity.eq("VALID_FOR_HYPOTHESIS"),"MODERATE","LOW")); lr["biological_reinterpretation"]=np.where(lr.ligand.eq("Pkm"),"Reclassify as intrinsic PKM-associated metabolic state; extracellular PKM→CD44 not established",np.where(lr.ligand.eq("Itgav"),"Integrin is represented in ligand position; direction is not supported",np.where(lr.ligand.isin(["Adam10","Adam17"]),"Membrane protease/shedding mechanism; specific substrate relationship required",""))); lr["valid_for_mechanistic_ranking"]=lr.interaction_validity.eq("VALID_FOR_HYPOTHESIS")
lr.to_csv(OUT/"interaction_biological_audit.tsv",sep="\t",index=False)
top=pd.concat([pd.read_csv(OUT/f"top_{c}.tsv",sep="\t").head(20) for c in ["extracellular","multicellular"]]).drop_duplicates("candidate").head(20); top=top.drop(columns=["interaction_validity","ligand_biological_class","receptor_biological_class","biological_reinterpretation","valid_for_mechanistic_ranking"],errors="ignore"); ta=top.merge(lr[["candidate","interaction_validity","ligand_biological_class","receptor_biological_class","biological_reinterpretation","valid_for_mechanistic_ranking"]],on="candidate",how="left"); ta["original_rank"]=range(1,len(ta)+1); ta["biological_validity"]=ta.interaction_validity; ta["ligand_validity"]=ta.ligand_biological_class; ta["receptor_validity"]=ta.receptor_biological_class; ta["direction_validity"]="ASSOCIATIVE_ONLY"; ta["mechanistic_plausibility"]=np.where(ta.valid_for_mechanistic_ranking,"SUPPORTED_AS_HYPOTHESIS","INSUFFICIENT"); ta["reclassification"]=ta.biological_reinterpretation; ta["reason"]=ta.interaction_validity; ta["recommended_action"]=np.where(ta.valid_for_mechanistic_ranking,"retain as hypothesis","exclude from LR ranking; retain appropriate intrinsic interpretation"); ta.to_csv(OUT/"top20_biological_audit.tsv",sep="\t",index=False)
(OUT/"interaction_biological_audit.md").write_text("# Interaction biological audit\n\nThis audit separates canonical/noncanonical extracellular signals, membrane-associated processing/contact proteins, metabolic enzymes and intracellular components. Presence in mouseconsensus/LIANA is not sufficient for mechanistic validity. PKM is classified as a metabolic intracellular protein and PKM→CD44 is excluded as an LR mechanism. ADAM10/17 are membrane proteases and are not described as soluble ligands. Unknown pairs remain unvalidated rather than promoted.\n")
print(lr.interaction_validity.value_counts().to_string())
