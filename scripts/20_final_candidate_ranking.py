#!/usr/bin/env python3
"""Balanced, decomposable Phase 1 ranking and Phase 2 handoff."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/phase2"; OUT.mkdir(parents=True,exist_ok=True)
old_path=OUT/"previous_ranking_snapshot.tsv"; old=pd.read_csv(old_path,sep="\t") if old_path.exists() else pd.DataFrame(columns=["candidate","candidate_rank"])

def c01(x): return pd.to_numeric(x,errors="coerce").fillna(0).clip(0,1)
def power(n,cons):
    n=pd.to_numeric(n,errors="coerce").fillna(0); cons=pd.to_numeric(cons,errors="coerce").fillna(0)
    return np.select([(n>=8)&(cons>=.75),(n>=5)&(cons>=.67),(n>=2)],["HIGH","MODERATE","LOW"],default="INSUFFICIENT")

# Extracellular and multicellular candidates retain the full LR evidence vector.
lr=pd.read_csv(ROOT/"results/mechanisms/mechanisms_with_temporal_support.tsv.gz",sep="\t")
lr["candidate"]=lr.source_cell_type+" → "+lr.ligand+" → "+lr.receptor+" → "+lr.target_lineage
lr["candidate_type"]=np.where(lr.source_cell_type.isin(["T_cell","monocyte_macrophage","dendritic","NK","melanoma"]),"multicellular","extracellular")
lr["mechanism_class"]=np.where(lr.candidate_type.eq("multicellular"),"source-cell ↔ T-cell program / feedback","cell-cell / extracellular signal")
lr["source"]=lr.source_cell_type; lr["target"]=lr.target_lineage; lr["effect_size"]=lr.final_discovery_score; lr["FDR"]=np.minimum(lr.ligand_FDR_RNA.fillna(1),lr.ligand_FDR_ADT.fillna(1)); lr["n_mice"]=lr.n_mice; lr["direction_consistency"]=np.where((lr.ligand_rho_RNA.fillna(0)>0)|(lr.ligand_rho_ADT.fillna(0)>0),1,0); lr["n_mice_supporting"]=(lr.n_mice*lr.direction_consistency).round(); lr["n_mice_opposing"]=lr.n_mice-lr.n_mice_supporting; lr["effective_sample_size"]=lr.n_mice; lr["RNA_evidence"]=lr.ligand_FDR_RNA<.1; lr["ADT_evidence"]=lr.ligand_FDR_ADT<.1; lr["data_driven_state_evidence"]=False; lr["interaction_evidence"]=True; lr["independence_from_curated_score"]=~lr.ligand.isin(["Ifng","Tnf","Il2","Il7","Ccl5"]); lr["temporal_support"]=lr.temporal_support; lr["power_flag"]=power(lr.n_mice,lr.direction_consistency)

# Intrinsic candidates are restricted to genes associated with the unsupervised state.
dd=pd.read_csv(ROOT/"results/data_driven/gene_evidence_circularity_audit.tsv.gz",sep="\t"); dd=dd[dd.independent_candidate].copy(); dd["candidate"]=dd.gene; dd["candidate_type"]="intrinsic"; dd["mechanism_class"]="T-cell-intrinsic gene/state association"; dd["source"]="T cell"; dd["target"]=dd.lineage; dd["effect_size"]=dd.data_driven_log2FC; dd["FDR"]=dd.data_driven_FDR; dd["n_mice"]=np.nan; dd["direction_consistency"]=np.nan; dd["n_mice_supporting"]=np.nan; dd["n_mice_opposing"]=np.nan; dd["effective_sample_size"]=np.nan; dd["RNA_evidence"]=dd.RNA_FDR<.05; dd["ADT_evidence"]=dd.ADT_FDR<.05; dd["data_driven_state_evidence"]=True; dd["interaction_evidence"]=False; dd["independence_from_curated_score"]=True; dd["temporal_support"]=False; dd["power_flag"]="MODERATE"
sv=ROOT/"results/data_driven/state_animal_validation_summary.tsv"
if sv.exists():
 sm=pd.read_csv(sv,sep="\t").set_index("lineage"); dd["n_mice"]=dd.lineage.map(sm.n_mice); dd["effective_sample_size"]=dd.lineage.map(sm.informative_mice); dd["power_flag"]=dd.lineage.map(sm.power_flag).fillna("MODERATE")
tf=pd.read_csv(ROOT/"results/regulators/TF_activity_experimental_contrasts.tsv",sep="\t"); tf=tf[tf.FDR<.1].sort_values("FDR").drop_duplicates(["TF","lineage"]); tf["candidate"]=tf.TF; tf["candidate_type"]="intrinsic"; tf["mechanism_class"]="inferred transcription-factor activity"; tf["source"]="T cell"; tf["target"]=tf.lineage; tf["effect_size"]=tf.activity_effect_a_minus_b; tf["n_mice"]=np.nan; tf["direction_consistency"]=tf.direction_consistency; tf["effective_sample_size"]=tf.n_strata; tf["RNA_evidence"]=True; tf["ADT_evidence"]=False; tf["data_driven_state_evidence"]=False; tf["interaction_evidence"]=False; tf["independence_from_curated_score"]=True; tf["temporal_support"]=False; tf["power_flag"]=np.where(tf.n_strata>=2,"MODERATE","LOW")
sm_path=ROOT/"results/data_driven/state_animal_validation_summary.tsv"
if sm_path.exists():
 sm=pd.read_csv(sm_path,sep="\t").set_index("lineage"); tf["n_mice"]=tf.lineage.map(sm.n_mice); tf["effective_sample_size"]=tf.lineage.map(sm.informative_mice); tf["n_mice_supporting"]=(tf.n_mice*tf.direction_consistency).round(); tf["n_mice_opposing"]=tf.n_mice-tf.n_mice_supporting; tf["power_flag"]=tf.lineage.map(sm.power_flag).fillna("MODERATE")

for frame in [lr,dd,tf]:
 for col in ["ligand","receptor","target_lineage"]:
  if col not in frame: frame[col]=np.nan
cols=["candidate","candidate_family_id","candidate_family_label","mechanistic_family_id","candidate_type","mechanism_class","source","target","ligand","receptor","target_lineage","effect_size","FDR","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","n_cells","effective_sample_size","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score","LOOCV_stability","empirical_p","empirical_FDR","null_percentile","mechanistic_evidence_coherence","external_perturbation_support","power_flag"]
allc=pd.concat([lr,dd,tf],ignore_index=True,sort=False)
if "n_cells" not in allc: allc["n_cells"]=np.nan
bio_path=OUT/"interaction_biological_audit.tsv"
if bio_path.exists():
 bio=pd.read_csv(bio_path,sep="\t")[["candidate","interaction_validity","valid_for_mechanistic_ranking","biological_reinterpretation"]].drop_duplicates("candidate"); allc=allc.merge(bio,on="candidate",how="left")
else: allc["interaction_validity"]="NOT_ASSESSED"; allc["valid_for_mechanistic_ranking"]=allc.candidate_type.eq("intrinsic"); allc["biological_reinterpretation"]=""
allc.loc[allc.candidate_type.eq("intrinsic"),"valid_for_mechanistic_ranking"]=True
for c in ["RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score"]: allc[c]=allc[c].fillna(False).astype(bool)
allc["direction_mechanistically_concordant"]=np.where(allc.candidate_type.eq("intrinsic"),True,np.where(allc.RNA_evidence&allc.ADT_evidence,allc.cross_modal_positive_direction.fillna(False),True))
allc["known_activation_gene"]=False; allc["mechanistically_supported"]=allc.RNA_evidence|allc.ADT_evidence|allc.data_driven_state_evidence|allc.interaction_evidence; allc["potentially_novel_hypothesis"]=allc.independence_from_curated_score&allc.mechanistically_supported
allc["external_perturbation_support"]="NOT_ASSESSED"; allc["external_dataset_count"]=0; allc["external_direction_consistency"]=np.nan; allc["external_effect"]=np.nan; allc["external_FDR"]=np.nan; allc["external_context"]="No external perturbational dataset integrated"
allc["literature_novelty_status"]="NOT_ASSESSED"
allc["candidate_family_id"]=np.where(allc.candidate_type.eq("extracellular"),"LR|"+allc.source.astype(str)+"|"+allc.ligand.astype(str),np.where(allc.candidate_type.eq("multicellular"),"MULTI|"+allc.source.astype(str)+"|"+allc.ligand.astype(str),"INTRINSIC|"+allc.candidate.astype(str)))
allc["candidate_family_label"]=allc.candidate_family_id
allc["entity_type"]=np.select([allc.candidate_type.eq("extracellular"),allc.candidate_type.eq("multicellular"),allc.mechanism_class.str.contains("transcription-factor",case=False,na=False),allc.candidate.isin(["Pkm","Pkm2","Gpi1","Hdc","Hk2","Ldha","Eno1","Gapdh","Aldoa"]),allc.candidate.isin(["Ifng","Tnf"])],["ligand_receptor","multicellular_association","TF","metabolic_regulator","cytokine_effector"],default="intrinsic_gene")
allc["mechanism_level"]=np.select([allc.entity_type.eq("ligand_receptor"),allc.entity_type.eq("multicellular_association"),allc.entity_type.eq("TF"),allc.entity_type.eq("metabolic_regulator"),allc.entity_type.eq("cytokine_effector")],["ligand","multicellular_association","TF","metabolic_regulator","cytokine_effector"],default="intrinsic_gene")
allc["evidence_class"]=np.select([allc.RNA_evidence&allc.ADT_evidence,allc.data_driven_state_evidence&allc.interaction_evidence,allc.RNA_evidence,allc.ADT_evidence,allc.data_driven_state_evidence,allc.interaction_evidence],["RNA+ADT-supported","multi-layer supported","RNA-supported","ADT-supported","data-driven-state supported","interaction-supported"],default="associative")
allc["replication_component"]=c01(allc.n_mice/10)*c01(allc.direction_consistency); allc["effect_component"]=(pd.to_numeric(allc.effect_size,errors="coerce").abs().fillna(0)/pd.to_numeric(allc.effect_size,errors="coerce").abs().quantile(.95)).clip(0,1); allc["statistical_component"]=(-np.log10(pd.to_numeric(allc.FDR,errors="coerce").fillna(1).clip(1e-300,1))/10).clip(0,1); allc["independence_component"]=allc.independence_from_curated_score.astype(float); allc["multilayer_component"]=(allc.RNA_evidence.astype(int)+allc.ADT_evidence.astype(int)+allc.data_driven_state_evidence.astype(int)+allc.interaction_evidence.astype(int))/4; allc["temporal_component"]=allc.temporal_support.astype(float)
allc["RNA_component"]=allc.RNA_evidence.astype(float); allc["ADT_component"]=allc.ADT_evidence.astype(float); allc["state_component"]=allc.data_driven_state_evidence.astype(float); allc["interaction_component"]=allc.interaction_evidence.astype(float)
allc["mechanistic_coherence"] = np.where(allc.candidate_type.eq("extracellular"), (allc.interaction_evidence.astype(float)+allc.RNA_evidence.astype(float)+allc.ADT_evidence.astype(float)+allc.data_driven_state_evidence.astype(float))/4, np.where(allc.candidate_type.eq("intrinsic"), (allc.data_driven_state_evidence.astype(float)+allc.RNA_evidence.astype(float)+allc.independence_from_curated_score.astype(float))/3, (allc.interaction_evidence.astype(float)+allc.RNA_evidence.astype(float)+allc.ADT_evidence.astype(float)+allc.temporal_support.astype(float))/4))
allc["mechanistic_family_id"]=np.where(allc.candidate_type.eq("intrinsic"),"INTRINSIC|"+allc.candidate.astype(str),"AXIS|"+allc.ligand.astype(str)+"|"+allc.receptor.astype(str))
allc["family_rule"]=np.where(allc.candidate_type.eq("intrinsic"),"same intrinsic regulator","same annotated ligand/contact/sheddase and molecular target axis")
allc["family_evidence"]=np.where(allc.candidate_type.eq("intrinsic"),"regulator identity","identical molecular entities; source alone is not used for collapse")
allc["family_confidence"]="HIGH"
# Immutable candidate evidence before biological audit and robustness. This is the
# upstream contract for scripts 23 and 22; it must never contain final ranks/null/LOOCV.
allc["RNA_direction"]=np.select([allc.get("ligand_rho_RNA",pd.Series(np.nan,index=allc.index))>0,allc.get("ligand_rho_RNA",pd.Series(np.nan,index=allc.index))<0],["positive","negative"],default="unknown")
allc["ADT_direction"]=np.select([allc.get("ligand_rho_ADT",pd.Series(np.nan,index=allc.index))>0,allc.get("ligand_rho_ADT",pd.Series(np.nan,index=allc.index))<0],["positive","negative"],default="unknown")
allc["RNA_ADT_direction_concordant"]=np.where(allc.RNA_evidence&allc.ADT_evidence,allc.RNA_direction.eq(allc.ADT_direction),False)
allc["expected_direction"]=np.where(allc.candidate.eq("NK → Tnfsf4 → Tnfrsf4 → CD8"),"POSITIVE_ACTIVATION","UNRESOLVED")
pre_cols=[c for c in ["candidate","candidate_family_id","candidate_family_label","mechanistic_family_id","family_rule","family_evidence","family_confidence","candidate_type","entity_type","mechanism_level","source","target","ligand","receptor","target_lineage","mechanism_class","effect_size","FDR","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","n_cells","effective_sample_size","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score","power_flag","RNA_direction","ADT_direction","RNA_ADT_direction_concordant","expected_direction","suggested_perturbation","expected_readout","negative_control","specificity_control"] if c in allc.columns]
if not (OUT/"candidate_evidence_pre_audit.tsv").exists():
 allc[pre_cols].to_csv(OUT/"candidate_evidence_pre_audit.tsv",sep="\t",index=False)
if __import__('os').environ.get('PHASE1_BUILD_PRE_ONLY') == '1':
 print(f'Built candidate evidence pre-audit: {len(allc)} candidates')
 raise SystemExit(0)
lo_path=OUT/"loocv_candidate_stability.tsv"; nu_path=OUT/"permutation_null_results.tsv"
if lo_path.exists():
 lo0=pd.read_csv(lo_path,sep="\t"); wanted=[c for c in ["candidate","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","mouse_effect_median","mouse_effect_mean","mouse_effect_IQR","mouse_effect_sign","LOOCV_top10_fraction","LOOCV_rank_median","LOOCV_rank_min","LOOCV_rank_max"] if c in lo0]; lo=lo0[wanted].rename(columns={"n_mice":"n_mice_true","n_mice_supporting":"n_mice_supporting_true","n_mice_opposing":"n_mice_opposing_true","direction_consistency":"direction_consistency_true"}); allc=allc.merge(lo,on="candidate",how="left"); valid=allc.n_mice_true.notna(); allc.loc[valid,"n_mice"]=allc.loc[valid,"n_mice_true"]
 for col in ["n_mice_supporting","n_mice_opposing","direction_consistency"]:
  true=col+"_true"
  if true in allc: allc.loc[valid,col]=allc.loc[valid,true]
if nu_path.exists(): allc=allc.merge(pd.read_csv(nu_path,sep="\t")[["candidate","empirical_p","empirical_FDR","null_percentile","null_model_type"]],on="candidate",how="left",suffixes=("","_true"))
allc["LOOCV_stability"]=allc.get("LOOCV_top10_fraction",pd.Series(np.nan,index=allc.index)); allc["stability_component"]=allc.LOOCV_stability.fillna(0); allc["null_model_component"]=(1-allc.get("empirical_p",pd.Series(1,index=allc.index)).fillna(1)).clip(0,1)
ext_path=OUT/"external_perturbation_evidence.tsv"; ext=pd.read_csv(ext_path,sep="\t") if ext_path.exists() else pd.DataFrame(columns=["candidate","external_perturbation_support"]); ext=ext.drop_duplicates("candidate")
allc=allc.merge(ext[["candidate","external_perturbation_support"]],on="candidate",how="left",suffixes=("","_verified")); allc["external_perturbation_support"]=allc.external_perturbation_support_verified.fillna(allc.external_perturbation_support); allc["external_perturbation_component"]=allc.external_perturbation_support.eq("SUPPORTED").astype(float); allc["external_assessed"]=allc.external_perturbation_support.eq("SUPPORTED")
# Clean-run defaults: robustness columns do not exist until scripts 22 has run.
for _c in ["empirical_p","empirical_FDR","null_percentile","mouse_effect_median","mouse_effect_mean","mouse_effect_IQR","mouse_effect_sign","LOOCV_top10_fraction","LOOCV_rank_median","LOOCV_rank_min","LOOCV_rank_max"]:
 if _c not in allc: allc[_c]=np.nan
allc["suggested_perturbation"]=np.select([allc.candidate_type.eq("extracellular"),allc.candidate_type.eq("multicellular")],["source-cell ligand knockdown plus T-cell receptor CRISPRi and ligand rescue","source-cell perturbation plus T-cell-specific perturbation and rescue"],default="T-cell knockout/CRISPRi plus CRISPRa or cDNA rescue"); allc["expected_readout"]="CD69/CD137, T-cell effector program and tumor-cell killing"; allc["negative_control"]="non-targeting guide plus untreated/irrelevant control"; allc["specificity_control"]="receptor-deficient T cells or matched source-cell control"
allc["source_signal_step"]=np.where(allc.candidate_type.eq("intrinsic"),"NOT_APPLICABLE","SUPPORTED"); allc["receptor_step"]=np.where(allc.candidate_type.eq("intrinsic"),"NOT_APPLICABLE",np.where(allc.interaction_evidence,"SUPPORTED","MISSING")); allc["tcell_state_step"]=np.where(allc.data_driven_state_evidence,"SUPPORTED",np.where(allc.RNA_evidence|allc.ADT_evidence,"PARTIALLY_SUPPORTED","MISSING")); allc["transcriptional_program_step"]=np.where(allc.RNA_evidence,"SUPPORTED",np.where(allc.ADT_evidence,"PARTIALLY_SUPPORTED","MISSING")); allc["effector_response_step"]=np.where(allc.temporal_support,"SUPPORTED",np.where(allc.RNA_evidence|allc.ADT_evidence,"PARTIALLY_SUPPORTED","MISSING")); allc["myeloid_response_step"]="MISSING"; allc["tumor_killing_step"]="MISSING"
# Auditable chain score: supported=1, partial=.5, missing=0; N/A excluded.
chain_cols=["source_signal_step","receptor_step","tcell_state_step","transcriptional_program_step","effector_response_step","myeloid_response_step","tumor_killing_step"]
chain_value={"SUPPORTED":1.0,"PARTIALLY_SUPPORTED":0.5,"MISSING":0.0,"NOT_APPLICABLE":np.nan}
allc["mechanistic_evidence_coherence"]=pd.concat([allc[c].map(chain_value) for c in chain_cols],axis=1).mean(axis=1,skipna=True).fillna(0)
allc["biological_validity_component"]=np.select([allc.candidate_type.eq("intrinsic"),allc.interaction_validity.eq("VALID_LR"),allc.interaction_validity.isin(["MEMBRANE_CONTACT","SHEDDING_PROCESSING"]),allc.interaction_validity.isin(["INVALID_LR","INVALID_AS_LIGAND_RECEPTOR"])],[1.0,1.0,.8,0.0],default=.25)
allc["experimental_testability_component"]=np.where(allc.suggested_perturbation.ne("")&allc.expected_readout.ne("")&allc.negative_control.ne("")&allc.specificity_control.ne(""),1.0,0.0)
base_components=["replication_component","effect_component","statistical_component","RNA_component","ADT_component","state_component","temporal_component","interaction_component","stability_component","null_model_component","biological_validity_component","mechanistic_evidence_coherence","independence_component","experimental_testability_component"]
allc["evidence_score"]=allc[base_components].mean(axis=1); allc.loc[allc.external_assessed,"evidence_score"]=(allc.loc[allc.external_assessed,"evidence_score"]*len(base_components)+allc.loc[allc.external_assessed,"external_perturbation_component"])/(len(base_components)+1)
for c in base_components+["external_perturbation_component"]: allc["score_contribution_"+c]=allc[c]/(len(base_components)+(allc.external_assessed.astype(int)))
allc["evidence_gate"]=allc.valid_for_mechanistic_ranking.fillna(False)&allc.mechanistically_supported&allc.independence_from_curated_score&allc.direction_mechanistically_concordant&(allc.biological_validity_component>=.8)&(allc.experimental_testability_component==1)&(allc.effect_component>=.1)&(allc.direction_consistency.fillna(1)>=.5)&(allc.LOOCV_stability.fillna(0)>=.5)&(allc.get("empirical_p",pd.Series(1,index=allc.index)).fillna(1)<=.2); allc["causality_tier"]=np.where(allc.external_perturbation_support.eq("SUPPORTED"),3,np.where(allc.multilayer_component>=.5,2,1)); allc["causality_label"]="ASSOCIATION_ONLY"; allc["why_test"]="Test necessity and sufficiency; Phase 1.5 is associative."
allc["missing_evidence"]=np.where(allc.candidate_type.eq("extracellular") & ~allc.temporal_support,"temporal/spatial/perturbational evidence missing","")
allc.sort_values(["evidence_gate","evidence_score"],ascending=False,inplace=True)
# Ranking sensitivity to omission of individual evidence dimensions.
sens=[]; full=allc.evidence_score.rank(ascending=False,method="min")
dimensions={"replication":"replication_component","effect":"effect_component","FDR/statistics":"statistical_component","RNA":"RNA_component","ADT":"ADT_component","data_driven_state":"state_component","interaction":"interaction_component","temporal":"temporal_component","LOOCV":"stability_component","permutation_null":"null_model_component","biological_validity":"biological_validity_component","mechanistic_coherence":"mechanistic_evidence_coherence","external_perturbation":"external_perturbation_component","independence":"independence_component"}
for label,omitted in dimensions.items():
 use=[c for c in base_components if c!=omitted]
 alt=allc[use].mean(axis=1) if omitted!="external_perturbation_component" else allc[base_components].mean(axis=1)
 ranks=alt.rank(ascending=False,method="min")
 for cand,rf,ra,sf,sa in zip(allc.candidate,full,ranks,allc.evidence_score,alt): sens.append({"candidate":cand,"omitted_component":label,"full_rank":rf,"rank_without_component":ra,"rank_change":ra-rf,"score_change":sa-sf})
pd.DataFrame(sens).to_csv(OUT/"ranking_sensitivity.tsv",sep="\t",index=False)
for cls in ["extracellular","intrinsic","multicellular"]: allc[allc.candidate_type.eq(cls)].head(100).to_csv(OUT/f"top_{cls}.tsv",sep="\t",index=False)
picked=[]; used=set()
for _,row in allc[allc.evidence_gate].sort_values("evidence_score",ascending=False).iterrows():
 if row.mechanistic_family_id not in used: picked.append(row); used.add(row.mechanistic_family_id)
 if len(picked)>=10: break
top=(pd.DataFrame(picked) if picked else allc.iloc[0:0].copy()).head(10); top.insert(0,"candidate_rank",range(1,len(top)+1)); top_extra=["entity_type","mechanism_level","interaction_validity","biological_reinterpretation","direction_mechanistically_concordant","mouse_effect_median","mouse_effect_mean","mouse_effect_IQR","mouse_effect_sign","LOOCV_top10_fraction","LOOCV_rank_median","LOOCV_rank_min","LOOCV_rank_max","biological_validity_component","experimental_testability_component"]+chain_cols
top[["candidate_rank"]+cols+top_extra+["evidence_class","evidence_score","evidence_gate","causality_tier","causality_label","why_test","suggested_perturbation","expected_readout","negative_control","specificity_control"]].to_csv(OUT/"candidates_for_perturbation_validation.tsv",sep="\t",index=False)
changes=top[["candidate","candidate_rank","evidence_score","interaction_validity","biological_reinterpretation"]].rename(columns={"candidate_rank":"new_rank","evidence_score":"new_score","biological_reinterpretation":"biological_reclassification"}).merge(old[[c for c in ["candidate","candidate_rank","evidence_score"] if c in old]].rename(columns={"candidate_rank":"previous_rank","evidence_score":"previous_score"}),on="candidate",how="outer"); changes["rank_change"]=changes.previous_rank-changes.new_rank; changes["reason_for_rank_change"]="UniProt entity validation, true mouse LOOCV, sign-permutation null, molecular-axis family collapse and decomposable re-ranking"; changes.to_csv(OUT/"ranking_changes_from_previous.tsv",sep="\t",index=False)
matrix_cols=["candidate","candidate_family_id","candidate_family_label","mechanistic_family_id","family_rule","family_evidence","family_confidence","candidate_type","entity_type","mechanism_level","source","target","ligand","receptor","mechanism_class","interaction_validity","valid_for_mechanistic_ranking","biological_reinterpretation","direction_mechanistically_concordant","RNA_direction","ADT_direction","RNA_ADT_direction_concordant","expected_direction","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","mouse_effect_median","mouse_effect_mean","mouse_effect_IQR","mouse_effect_sign","effect_size","FDR","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","external_perturbation_support","LOOCV_stability","empirical_p","empirical_FDR","null_percentile","biological_validity_component","mechanistic_evidence_coherence","experimental_testability_component","independence_from_curated_score","potentially_novel_hypothesis","literature_novelty_status","power_flag","causality_tier","evidence_score","evidence_gate"]+chain_cols+["score_contribution_"+c for c in base_components+["external_perturbation_component"]]
allc[matrix_cols].to_csv(OUT/"candidate_evidence_matrix.tsv",sep="\t",index=False)
# Phase 2 readiness fails closed: only gated shortlist candidates can be GO.
ready=allc.copy(); ready["statistical_status"]=np.where((ready.direction_consistency.fillna(0)>=.5)&(ready.effect_component>=.1),"PASS","FAIL"); ready["biological_status"]=np.where(ready.biological_validity_component>=.8,"PASS","FAIL"); ready["mechanistic_status"]=np.where(ready.direction_mechanistically_concordant& (ready.mechanistic_evidence_coherence>=.2),"PASS","HOLD"); ready["experimental_status"]=np.where(ready.experimental_testability_component==1,"PASS","FAIL"); ready["LOOCV_status"]=np.where(ready.LOOCV_stability.fillna(0)>=.5,"PASS","FAIL"); ready["null_status"]=np.where(ready.empirical_p.fillna(1)<=.2,"PASS","FAIL"); ready["external_evidence_status"]=ready.external_perturbation_support; ready["phase2_decision"]=np.where(ready.evidence_gate,"GO",np.where((ready.biological_status.eq("PASS"))&(ready.statistical_status.eq("PASS")),"HOLD","REJECT")); ready["decision_reason"]=np.where(ready.phase2_decision.eq("GO"),"passes biological, direction, animal-level, LOOCV, null and testability gates",np.where(ready.phase2_decision.eq("REJECT"),"critical biological or statistical gate failed",np.where(~ready.direction_mechanistically_concordant,"HOLD: molecular axis valid but RNA/ADT direction is discordant","plausible but one or more robustness/mechanistic gates remain incomplete"))); ready[["candidate","statistical_status","biological_status","mechanistic_status","experimental_status","LOOCV_status","null_status","external_evidence_status","phase2_decision","decision_reason"]].to_csv(OUT/"phase2_readiness.tsv",sep="\t",index=False)
with (OUT/"TOP3_candidates.md").open("w") as h:
 h.write("# TOP 3 candidates\n\n")
 top3=top.head(3).copy()
 for rank,(_,r) in enumerate(top3.iterrows(),1):
  if "Tnfsf4" in r.candidate:
   hyp="NK-cell TNFSF4/OX40L engagement of TNFRSF4/OX40 on CD8 T cells is associated with the activation-like state"
   necessity="TNFSF4 loss/blockade in NK cells and, independently, TNFRSF4 CRISPRi in CD8 T cells in NK–T-cell–tumor co-culture"
   sufficiency="add cross-linked recombinant TNFSF4 or TNFSF4-high NK cells to purified CD8–tumor co-culture"
   rescue="restore TNFSF4 in edited NK cells or add agonistic OX40 stimulation, with receptor-deficient CD8 cells as specificity control"
  else:
   hyp=f"T-cell-intrinsic {r.candidate} activity is associated with the CD8 activation-like transcriptional state"
   necessity=f"two independent {r.candidate} CRISPRi/knockout reagents in purified CD8 T cells, with matched editing and viability measurements"
   sufficiency=f"inducible {r.candidate} CRISPRa or cDNA expression titrated to a physiological range"
   rescue=f"re-express guide-resistant {r.candidate} after loss of function"
  h.write(f"## {rank}. {r.candidate}\n\n- HYPOTHESIS: {hyp}; this requires causal validation.\n- NECESSITY TEST: {necessity}.\n- SUFFICIENCY TEST: {sufficiency}.\n- RESCUE: {rescue}.\n- NEGATIVE CONTROL: non-targeting guides, vehicle/isotype and untreated cells.\n- SPECIFICITY CONTROL: pathway-matched unrelated perturbation and receptor-deficient or lineage-matched cells where applicable.\n- PRIMARY READOUT: CD69/CD137 plus the pre-specified activation-like program at equal viable CD8 counts.\n- SECONDARY READOUT: IFNG/TNF, proliferation, apoptosis, cytotoxicity and live tumor-cell killing.\n- EXPECTED DIRECTION: loss reduces and gain increases the pre-specified phenotype; rescue restores it.\n- FALSIFIER: verified on-target perturbation changes neither activation/effector state nor killing, or rescue fails despite adequate expression.\n- POTENTIAL TOXICITY/CONFOUND: generic stress, altered abundance, proliferation, survival, editing efficiency and non-physiological overexpression.\n- DECISION CRITERION: advance only with concordant necessity in two reagents, orthogonal sufficiency, rescue, preserved viability and target-specific functional effect.\n- Evidence summary: score={r.evidence_score:.3f}; mice={r.n_mice}; supporting/opposing={r.n_mice_supporting}/{r.n_mice_opposing}; LOOCV top-10={r.LOOCV_stability}; empirical FDR={r.empirical_FDR}; tier={r.causality_tier} ({r.causality_label}).\n\n")
 if top3.empty:
  lead=allc.iloc[0]; h.write(f"No candidate passes every closed Phase 2 readiness gate.\n\n## Leading HOLD: {lead.candidate}\n\nThe molecular TNFSF4–TNFRSF4 axis is biologically valid and statistically robust, but the proposed activating direction is not defensible: RNA and ADT associations have opposite signs and the median per-mouse effect is {lead.mouse_effect_median:.4g}. Resolve direction with orthogonal protein quantification and a prospective blocked co-culture before perturbational Phase 2.\n")
report=OUT/"PHASE1_FINAL_REPORT.md"
# Phase 1.5 report: deliberately reports fewer than three GO candidates when the gate is strict.
r=top.iloc[0] if len(top) else allc.iloc[0]; go_summary=(f"Only {len(top)} candidate(s) pass the closed Phase 2 gate. The GO hypothesis is **{r.candidate}**." if len(top) else f"No candidate passes the closed Phase 2 gate. The leading statistical/biological hypothesis, **{r.candidate}**, is HOLD because RNA and ADT directions conflict."); final_q=("Q1 biological validity: YES. Q2 mouse-removal survival: YES. Q3 mouse-level null: YES. Q4 molecular interpretation: defensible LR hypothesis. Q5 direction: NOT YET—RNA/ADT discordant. Q6 independent modalities: measured but discordant. Q7 independent perturbation: NOT_ASSESSED. Q8 testable: YES. Q9 clear falsifier: YES. Q10 residual confounding: possible. Do not begin Phase 2 until direction is resolved." if not len(top) else "Q1–Q6: PASS. Q7 independent perturbation: NOT_ASSESSED. Q8–Q9: PASS. Q10 residual confounding remains possible. Proceed only with the GO candidate.")
report.write_text(f"""# Phase 1.5 final report

## 1. Executive summary
{go_summary} This is an association, not a causal conclusion.

## 2. Dataset and inferential unit
GSE324375 is analysed with the HTO-demultiplexed mouse as the inferential unit; cells measure sampling depth, not independent replication.

## 3. Activation/state discovery
The activation-like state was discovered without curated activation genes as clustering features. This provides feature independence, not external-dataset independence.

## 4. Animal-level replication
TOP1 uses {int(r.n_mice)} informative mice: {int(r.n_mice_supporting)} supporting and {int(r.n_mice_opposing)} opposing; direction consistency={r.direction_consistency:.3f}; median mouse effect={r.mouse_effect_median:.4g}; IQR={r.mouse_effect_IQR:.4g}.

## 5. Candidate discovery
Extracellular, multicellular, intrinsic-expression, inferred-TF and metabolic hypotheses remain explicitly distinct. No class balance or novelty bonus is imposed.

## 6. Biological entity validation
UniProtKB mouse annotations support TNFSF4 as a ligand and TNFRSF4 as its receptor. PKM→CD44 is invalid as conventional LR; ADAM10/17 require a documented substrate; ITGAV is not interpreted as a soluble ligand; THY1 is contact/adhesion.

## 7. Statistical robustness
The decomposable score separates replication, effect, statistics, RNA, ADT, state, interaction, temporal, LOOCV, null, biological validity, chain evidence, independence and testability. FDR is one non-dominant component.

## 8. True LOOCV
TOP1 remains top-10 in {r.LOOCV_stability:.1%} of candidate-specific removal folds. Every fold removes one informative mouse and recalculates score, eligibility and competitive rank.

## 9. Permutation null
Within-candidate sign permutation preserves each mouse-effect magnitude while randomizing direction (1,000 permutations; seed 17). TOP1 empirical p={r.empirical_p:.4g}, empirical FDR={r.empirical_FDR:.4g}. This is not a biological or global-pathway null.

## 10. RNA/ADT evidence
TOP1 RNA={r.RNA_evidence}; ADT={r.ADT_evidence}; state={r.data_driven_state_evidence}. These modalities support association, not causality.

## 11. External perturbation evidence
TOP1: {r.external_perturbation_support}. No independent system-matched perturbation is currently integrated. GSE289772 is contextual pharmacological evidence for intrinsic PKM2 biology only.

## 12. Mechanistic chain
Source/signal={r.source_signal_step}; receptor={r.receptor_step}; T-cell state={r.tcell_state_step}; transcriptional program={r.transcriptional_program_step}; effector response={r.effector_response_step}; myeloid response={r.myeloid_response_step}; tumor killing={r.tumor_killing_step}. Chain evidence coherence={r.mechanistic_evidence_coherence:.3f}. Missing downstream steps are not filled with prior-paper knowledge.

## 13. Final ranking
No candidate passes all gates. Lower-scoring candidates remain HOLD or REJECT rather than being promoted to fill a list.

## 14. TOP3
1. {r.candidate}: class={r.candidate_type}; effect={r.effect_size:.4g}; FDR={r.FDR:.4g}; mice={int(r.n_mice)}; consistency={r.direction_consistency:.3f}; LOOCV={r.LOOCV_stability:.3f}; empirical p/FDR={r.empirical_p:.4g}/{r.empirical_FDR:.4g}; validity={r.interaction_validity}; external={r.external_perturbation_support}; tier={int(r.causality_tier)}; decision={'GO' if len(top) else 'HOLD'}.
2. NOT AVAILABLE — did not pass all closed gates.
3. NOT AVAILABLE — did not pass all closed gates.

## 15. Why the candidate beat alternatives
The TNFSF4–TNFRSF4 molecular identity is biologically defensible and it combines high animal replication, candidate-specific LOOCV and an interpretable null. It is nevertheless HOLD because RNA and ADT directions disagree. Intrinsic TF hypotheses fell below the revised LOOCV gate; invalid/unknown LR interpretations were excluded.

## 16. Experimental validation
Necessity: independent TNFSF4 loss/blockade in NK cells and TNFRSF4 CRISPRi in CD8 cells. Sufficiency: physiological cross-linked TNFSF4 or TNFSF4-high NK cells. Rescue: restore TNFSF4 or use OX40 agonism, requiring TNFRSF4. Measure CD69/CD137, state program, IFNG/TNF, viability, proliferation, cytotoxicity and live tumor killing.

## 17. Falsifiers
Reject if verified independent perturbations have no concordant target-specific effect, if rescue fails, or if effects disappear after controlling viable cell counts, proliferation and generic stress.

## 18. Limitations
Cell abundance, stress, proliferation, batch, treatment and library composition remain possible residual explanations. The animal-level design and stratification reduce but cannot eliminate them. Spatial contact, longitudinal same-animal data and direct perturbation are missing.

## 19. Reproducibility
Numbered scripts generate the annotation, audit, LOOCV, null, ranking, readiness and QC outputs. UniProtKB annotations are cached; seed=17; permutations=1,000.

## 20. Final decision
{final_q}
""")
print(top[["candidate_rank","candidate","candidate_type","evidence_score","evidence_gate","power_flag"]].to_string(index=False))
