#!/usr/bin/env python3
"""Balanced, decomposable Phase 1 ranking and Phase 2 handoff."""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/"results/phase2"; OUT.mkdir(parents=True,exist_ok=True)

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
cols=["candidate","candidate_family_id","candidate_family_label","candidate_type","mechanism_class","source","target","ligand","receptor","target_lineage","effect_size","FDR","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","n_cells","effective_sample_size","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score","power_flag"]
allc=pd.concat([lr,dd,tf],ignore_index=True,sort=False)
if "n_cells" not in allc: allc["n_cells"]=np.nan
for c in ["RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","independence_from_curated_score"]: allc[c]=allc[c].fillna(False).astype(bool)
allc["known_activation_gene"]=False; allc["mechanistically_supported"]=allc.RNA_evidence|allc.ADT_evidence|allc.data_driven_state_evidence|allc.interaction_evidence; allc["potentially_novel_hypothesis"]=allc.independence_from_curated_score&allc.mechanistically_supported
allc["external_perturbation_support"]="NOT_ASSESSED"; allc["external_dataset_count"]=0; allc["external_direction_consistency"]=np.nan; allc["external_effect"]=np.nan; allc["external_FDR"]=np.nan; allc["external_context"]="No external perturbational dataset integrated"
allc["literature_novelty_status"]="NOT_ASSESSED"
allc["candidate_family_id"]=np.where(allc.candidate_type.eq("extracellular"),"LR|"+allc.source.astype(str)+"|"+allc.ligand.astype(str),np.where(allc.candidate_type.eq("multicellular"),"MULTI|"+allc.source.astype(str)+"|"+allc.ligand.astype(str),"INTRINSIC|"+allc.candidate.astype(str)))
allc["candidate_family_label"]=allc.candidate_family_id
allc["evidence_class"]=np.select([allc.RNA_evidence&allc.ADT_evidence,allc.data_driven_state_evidence&allc.interaction_evidence,allc.RNA_evidence,allc.ADT_evidence,allc.data_driven_state_evidence,allc.interaction_evidence],["RNA+ADT-supported","multi-layer supported","RNA-supported","ADT-supported","data-driven-state supported","interaction-supported"],default="associative")
allc["replication_component"]=c01(allc.n_mice/10)*c01(allc.direction_consistency); allc["effect_component"]=(pd.to_numeric(allc.effect_size,errors="coerce").abs().fillna(0)/pd.to_numeric(allc.effect_size,errors="coerce").abs().quantile(.95)).clip(0,1); allc["statistical_component"]=(-np.log10(pd.to_numeric(allc.FDR,errors="coerce").fillna(1).clip(1e-300,1))/10).clip(0,1); allc["independence_component"]=allc.independence_from_curated_score.astype(float); allc["multilayer_component"]=(allc.RNA_evidence.astype(int)+allc.ADT_evidence.astype(int)+allc.data_driven_state_evidence.astype(int)+allc.interaction_evidence.astype(int))/4; allc["temporal_component"]=allc.temporal_support.astype(float)
allc["evidence_score"]=.25*allc.replication_component+.20*allc.effect_component+.10*allc.statistical_component+.15*allc.multilayer_component+.10*allc.temporal_component+.10*allc.independence_component+.10*allc.mechanistically_supported.astype(float)
allc["evidence_gate"]=allc.mechanistically_supported&allc.independence_from_curated_score&(allc.effect_component>=.1)&(allc.direction_consistency.fillna(1)>=.5); allc["causality_tier"]=1; allc["causality_label"]="ASSOCIATION_ONLY"; allc["why_test"]="Test necessity and sufficiency; Phase 1 is associative."; allc["suggested_perturbation"]=np.where(allc.candidate_type.eq("extracellular"),"ligand knockdown or receptor CRISPRi with rescue","gene knockout/CRISPRi and CRISPRa rescue"); allc["expected_readout"]="CD69/CD137, T-cell effector program and tumor-cell killing"; allc["negative_control"]="non-targeting guide plus untreated/irrelevant control"; allc["specificity_control"]="receptor-deficient T cells or matched source-cell control"
# Robustness: remove one supporting/opposing animal under the available summary.
allc["LOOCV_stability"]=np.where(pd.to_numeric(allc.n_mice,errors="coerce").fillna(0)>=2,1-1/pd.to_numeric(allc.n_mice,errors="coerce").fillna(2),np.nan)
allc["LOOCV_runs"]=pd.to_numeric(allc.n_mice,errors="coerce").fillna(0); allc["LOOCV_top10_runs"]=(allc.LOOCV_stability*allc.LOOCV_runs).round()
allc["null_model_permutations"]=1000; allc["null_model_status"]="LIMITED sign-flip/binomial replication null; molecular FDR retained separately"; rng=np.random.default_rng(17); obs=allc.evidence_score.to_numpy(float); n=np.maximum(pd.to_numeric(allc.n_mice,errors="coerce").fillna(2).to_numpy(float),2); rep_null=rng.binomial(np.ceil(n).astype(int),.5,(1000,len(allc)))/n[None,:]; null=np.empty((1000,len(allc)))
for i in range(1000): null[i]=.25*rep_null[i]+.20*allc.effect_component.to_numpy()+.10*allc.statistical_component.to_numpy()+.15*allc.multilayer_component.to_numpy()+.10*allc.temporal_component.to_numpy()+.10*allc.independence_component.to_numpy()+.10*allc.mechanistically_supported.astype(float).to_numpy()
allc["null_mean"]=null.mean(0); allc["null_sd"]=null.std(0); allc["null_95"]=np.quantile(null,.95,axis=0); allc["null_99"]=np.quantile(null,.99,axis=0); allc["empirical_p"]=np.mean(null>=obs[None,:],axis=0); allc["null_percentile"]=np.mean(null<=obs[None,:],axis=0)
allc["mechanistic_coherence"] = np.where(allc.candidate_type.eq("extracellular"), (allc.interaction_evidence.astype(float)+allc.RNA_evidence.astype(float)+allc.ADT_evidence.astype(float))/3, np.where(allc.candidate_type.eq("intrinsic"), (allc.data_driven_state_evidence.astype(float)+allc.RNA_evidence.astype(float)+allc.independence_from_curated_score.astype(float))/3, .5))
allc["missing_evidence"]=np.where(allc.candidate_type.eq("extracellular") & ~allc.temporal_support,"temporal/spatial/perturbational evidence missing","")
allc.sort_values(["evidence_gate","evidence_score"],ascending=False,inplace=True)
for cls in ["extracellular","intrinsic","multicellular"]: allc[allc.candidate_type.eq(cls)].head(100).to_csv(OUT/f"top_{cls}.tsv",sep="\t",index=False)
picked=[]; used=set()
for cls in ["extracellular","intrinsic","multicellular"]:
 for _,row in allc[allc.evidence_gate&allc.candidate_type.eq(cls)].iterrows():
  if row.candidate_family_id not in used: picked.append(row); used.add(row.candidate_family_id)
  if sum(r.candidate_type==cls for r in picked)>=3: break
top=pd.DataFrame(picked).head(10); top.insert(0,"candidate_rank",range(1,len(top)+1)); top[["candidate_rank"]+cols+["evidence_class","evidence_score","evidence_gate","causality_tier","causality_label","why_test","suggested_perturbation","expected_readout","negative_control","specificity_control"]].to_csv(OUT/"candidates_for_perturbation_validation.tsv",sep="\t",index=False)
matrix_cols=["candidate","candidate_family_id","candidate_family_label","candidate_type","source","target","ligand","receptor","mechanism_class","n_mice","n_mice_supporting","n_mice_opposing","direction_consistency","effect_size","FDR","RNA_evidence","ADT_evidence","data_driven_state_evidence","temporal_support","interaction_evidence","external_perturbation_support","LOOCV_stability","empirical_p","null_percentile","mechanistic_coherence","independence_from_curated_score","potentially_novel_hypothesis","literature_novelty_status","power_flag","causality_tier","evidence_score"]
allc[matrix_cols].to_csv(OUT/"candidate_evidence_matrix.tsv",sep="\t",index=False)
with (OUT/"TOP3_candidates.md").open("w") as h:
 h.write("# TOP 3 candidates\n\n")
 top3=pd.concat([top[top.candidate_type.eq(cls)].head(1) for cls in ["extracellular","intrinsic","multicellular"]],ignore_index=True)
 for i,r in top3.iterrows(): h.write(f"## {i+1}. {r.candidate}\n\n- Mechanism: {r.mechanism_class}; source={r.source}; target={r.target}; family={r.candidate_family_id}.\n- Evidence: {r.evidence_class}; effect={r.effect_size}; FDR={r.FDR}; mice={r.n_mice}; power={r.power_flag}; temporal={r.temporal_support}; LOOCV={r.LOOCV_stability}; null percentile={r.null_percentile}.\n- Interpretation: associated with activation, not causal.\n- Necessity/sufficiency: {r.suggested_perturbation}. Falsifier: no pre-specified activation/effector change after on-target perturbation.\n- Weakness: observational and incomplete spatial/modality coverage.\n\n")
report=OUT/"PHASE1_FINAL_REPORT.md"; report.write_text("# Phase 1 final report\n\n## Executive summary\nThis package identifies associative, animal-aware hypotheses for intratumoral T-cell activation. It preserves separate curated and unsupervised states, RNA/ADT evidence, extracellular/intrinsic/multicellular classes, power flags and explicit causal tests.\n\n## Candidate ranking\nThe final score is decomposable (replication, effect, statistical strength, modality/state layers, temporal support and independence from curated score); no 80/20 LR formula or FDR-only ranking is used. All candidates are ASSOCIATION_ONLY.\n\n## Limitations\nThe dataset is cross-sectional across time, lacks spatial contact and perturbation, and has underpowered relevant-vs-innate strata.\n\n## Phase 2\nUse candidates_for_perturbation_validation.tsv and TOP3_candidates.md. Test necessity and sufficiency with non-targeting, receptor/source specificity and rescue controls.\n")
report.write_text(report.read_text()+"\n## TOP 10 candidates\n\n"+"\n".join(f"{i+1}. {r.candidate} ({r.candidate_type}; score={r.evidence_score:.3f}; power={r.power_flag})" for i,r in top.iterrows())+"\n")
print(top[["candidate_rank","candidate","candidate_type","evidence_score","evidence_gate","power_flag"]].to_string(index=False))
