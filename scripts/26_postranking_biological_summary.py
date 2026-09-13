#!/usr/bin/env python3
"""Create presentation-only top-20 biological audit after final ranking."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]; out=root/'results/phase2'; src=out/'candidate_evidence_matrix.tsv'; bio=out/'interaction_biological_audit.tsv'
if src.exists() and bio.exists():
 frame=pd.read_csv(src,sep='\t',low_memory=False); score_col='final_score' if 'final_score' in frame else 'evidence_score'; top=frame.sort_values(score_col,ascending=False).head(20)[['candidate']]
 b=pd.read_csv(bio,sep='\t').merge(top,on='candidate',how='inner'); b.to_csv(out/'top20_biological_audit.tsv',sep='\t',index=False)
 eff=out/'candidate_mouse_effects.tsv'; lo=out/'loocv_candidate_stability.tsv'
 if eff.exists() and lo.exists():
  e=pd.read_csv(eff,sep='\t'); s=e.groupby('candidate').effect_value.agg(n_positive=lambda x:int((x>0).sum()),n_negative=lambda x:int((x<0).sum()),n_zero_or_undetermined=lambda x:int((x==0).sum()),n_mice_effects='count').reset_index(); s['majority_direction']=s.apply(lambda r:'positive' if r.n_positive>r.n_negative else ('negative' if r.n_negative>r.n_positive else 'balanced'),axis=1); s['majority_direction_fraction']=s[['n_positive','n_negative']].max(axis=1)/s.n_mice_effects; s['two_sided_direction_strength']=2*(s.majority_direction_fraction-.5).abs(); l=pd.read_csv(lo,sep='\t'); l=l.drop(columns=[c for c in ['n_positive','n_negative','n_zero_or_undetermined','majority_direction','majority_direction_fraction','two_sided_direction_strength','n_mice_effects'] if c in l.columns]); z=l.merge(s,on='candidate',how='left'); z=z.drop(columns=[c for c in ['n_mice_x','n_mice_y'] if c in z.columns]).rename(columns={'n_mice_x':'n_mice'}); z.to_csv(lo,sep='\t',index=False)
 print(f'Post-ranking biological summary: {len(b)} candidates')
# Normalize merge suffixes so the LOOCV table retains a single n_mice column.
if lo.exists():
 l=pd.read_csv(lo,sep='\t')
 if 'n_mice' not in l.columns and {'n_positive','n_negative','n_zero_or_undetermined'}.issubset(l.columns):
  l['n_mice']=l['n_positive']+l['n_negative']+l['n_zero_or_undetermined']
 if 'n_mice' not in l.columns:
  for c in ('n_mice_x','n_mice_y'):
   if c in l.columns:
    l['n_mice']=l[c]; break
  l=l.drop(columns=[c for c in ('n_mice_x','n_mice_y') if c in l.columns])
  l.to_csv(lo,sep='\t',index=False)
 else:
  l.to_csv(lo,sep='\t',index=False)
# Synchronize explicit leading-candidate headers from readiness + ranking.
ready_path=out/'phase2_readiness.tsv'
if src.exists() and ready_path.exists():
 ranking=pd.read_csv(src,sep='\t',low_memory=False); ready=pd.read_csv(ready_path,sep='\t');
 def lead(dec):
  c=ready.loc[ready.phase2_decision.eq(dec),'candidate']; q=ranking[ranking.candidate.isin(c)].sort_values(['evidence_score','candidate'],ascending=[False,True]); return 'NONE' if q.empty else str(q.iloc[0].candidate)
 lg,lh=lead('GO'),lead('HOLD'); report=out/'PHASE1_FINAL_REPORT.md'; old=report.read_text() if report.exists() else ''; report.write_text(f'Leading GO: {lg}\nLeading HOLD: {lh}\n\n'+old)
 holds=ready[ready.phase2_decision.eq('HOLD')][['candidate']].merge(ranking,on='candidate').sort_values(['evidence_score','candidate'],ascending=[False,True]).head(3)
 lines=['# Phase 1.5 final report','',f'Leading GO: {lg}',f'Leading HOLD: {lh}','',('No candidate passes every closed Phase 2 readiness gate.' if lg=='NONE' else f'Leading GO candidate: {lg}'),'','## Leading HOLD hypotheses']
 for i,(_,row) in enumerate(holds.iterrows(),1): lines.append(f"{i}. {row.candidate} — class={row.candidate_type}; score={row.evidence_score:.4f}; decision=HOLD.")
 for i in range(len(holds)+1,4): lines.append(f'{i}. NOT AVAILABLE — no additional candidate selected.')
 top=ranking[ranking.candidate.eq(lh)].iloc[0] if lh!='NONE' else ranking.iloc[0]; top_lo=pd.read_csv(lo,sep='\t'); top_lo=top_lo[top_lo.candidate.eq(lh)].iloc[0] if lh!='NONE' and (top_lo.candidate==lh).any() else None
 lines += ['', '## 1. Executive summary','',f'No candidate passes every closed Phase 2 readiness gate. The leading HOLD hypothesis is **{lh}**, classified as {top.candidate_type}. This is an observational prioritization and does not establish causality or therapeutic benefit.','', '## 2. Dataset and inferential unit','', 'GSE324375 is analysed with the HTO-demultiplexed mouse as the inferential unit; cells represent sampling depth, not independent replication.','', '## 3. Leading candidate','',f'Candidate: {lh}',f'Candidate class: {top.candidate_type}',f'Evidence score: {top.evidence_score:.4f}','Phase 2 decision: HOLD','']
 if top_lo is not None:
  lines += ['## 4. Animal-level replication','',f"Informative mice: {int(top_lo.n_mice)}.",f"Expected causal direction: {top.expected_direction}.","Hypothesis-supporting mice: NOT APPLICABLE." if top.expected_direction=='UNRESOLVED' else f"Hypothesis-supporting mice: {int(top_lo.n_mice_supporting)}.","Hypothesis-opposing mice: NOT APPLICABLE." if top.expected_direction=='UNRESOLVED' else f"Hypothesis-opposing mice: {int(top_lo.n_mice_opposing)}.",f"Positive mouse effects: {int(top_lo.n_positive)}.",f"Negative mouse effects: {int(top_lo.n_negative)}.",f"Zero/undetermined effects: {int(top_lo.n_zero_or_undetermined)}.",f"Majority direction: {top_lo.majority_direction}.",f"Majority direction fraction: {top_lo.majority_direction_fraction:.3f}.",f"Two-sided direction strength: {top_lo.two_sided_direction_strength:.3f}.",f"Median mouse effect: {top.mouse_effect_median:.4f}.",f"IQR: {top.mouse_effect_IQR:.4f}.",'', '## 5. LOOCV robustness','',f'LOOCV stability: {top.LOOCV_stability:.3f}.','Each candidate-specific fold removes one informative mouse and recalculates eligibility and competitive ranking.','', '## 6. Permutation null','',f"Empirical p-value: {top.empirical_p}.",f"Empirical FDR: {top.empirical_FDR}.",f"Null percentile: {top.null_percentile}.",'The null is candidate-specific and does not constitute a global-pathway or causal null.','']
 lines += ['## 7. Multimodal evidence','',f"RNA evidence: {top.RNA_evidence}; ADT evidence: {top.ADT_evidence}; state evidence: {top.data_driven_state_evidence}.",'', 'The association remains observational; independent perturbation is required before Phase 2.']
 lines += ['## 8. External perturbation evidence','',f"Status: {getattr(top,'external_perturbation_support','NOT AVAILABLE') or 'NOT AVAILABLE'}.",'Independent system-matched perturbation evidence is required before causal interpretation.','', '## 9. Biological validity','',f"Entity type: {getattr(top,'entity_type','NOT AVAILABLE')}.",f"Mechanism class: {getattr(top,'mechanism_class','NOT AVAILABLE')}.",f"Interaction validity: {getattr(top,'interaction_validity','NOT AVAILABLE')}.",'', '## 10. Mechanistic chain','']
 for c in [c for c in ranking.columns if c.startswith('chain_')]: lines.append(f"{c.replace('chain_','').replace('_',' ').title()}: {top[c]}.")
 lines += ['', '## 11. Leading HOLD hypotheses',''] + [f"{i}. {row.candidate} — class={row.candidate_type}; score={row.evidence_score:.4f}; decision=HOLD." for i,(_,row) in enumerate(holds.iterrows(),1)] + ['', '## 12. Why the leading candidate remains HOLD','', f"Decision reason: {ready.loc[ready.candidate.eq(lh),'decision_reason'].iloc[0] if (ready.candidate==lh).any() else 'NOT AVAILABLE'}.",'', '## 13. Experimental validation','', 'Use candidate-class-specific perturbation, matched controls, rescue where feasible, and measure activation, viability, abundance and tumor-cell killing.','', '## 14. Falsifiers','', '- No reproducible phenotype after independent perturbation.','- Rescue fails.','- Effect is explained by viability or abundance.','', '## 15. Limitations','', 'The analysis is associative; residual confounding and lack of system-matched perturbation remain.']
 else:
  lines += ['## 4. Animal-level replication','', 'Informative-mouse summary: NOT AVAILABLE.','Reason: the leading candidate is not represented in the predeclared LOOCV/mouse-effect summary universe.','', '## 5. LOOCV robustness','', 'LOOCV stability: NOT AVAILABLE.','Reason: candidate-specific folds are unavailable.','', '## 6. Permutation null','', 'Empirical p-value: NOT AVAILABLE.','Empirical FDR: NOT AVAILABLE.','Null percentile: NOT AVAILABLE.','The null is candidate-specific and does not establish biological causality.','', '## 7. Multimodal evidence','',f"RNA evidence: {top.RNA_evidence}.",f"ADT evidence: {top.ADT_evidence}.",f"Data-driven state evidence: {top.data_driven_state_evidence}.",'These modalities support prioritization but do not establish causality.']
 report.write_text('\n'.join(lines)+'\n')
 (out/'TOP3_candidates.md').write_text('\n'.join(['# TOP 3 candidates','',f'Leading GO: {lg}',f'Leading HOLD: {lh}','']+[f"{i}. {row.candidate} — class={row.candidate_type}; score={row.evidence_score:.4f}; decision=HOLD." for i,(_,row) in enumerate(holds.iterrows(),1)]+[f'{i}. NOT AVAILABLE — no additional candidate selected.' for i in range(len(holds)+1,4)])+'\n')
