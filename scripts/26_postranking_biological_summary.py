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
