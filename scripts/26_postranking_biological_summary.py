#!/usr/bin/env python3
"""Create presentation-only top-20 biological audit after final ranking."""
from pathlib import Path
import pandas as pd
root=Path(__file__).resolve().parents[1]; out=root/'results/phase2'; src=out/'candidate_evidence_matrix.tsv'; bio=out/'interaction_biological_audit.tsv'
if src.exists() and bio.exists():
 frame=pd.read_csv(src,sep='\t',low_memory=False); score_col='final_score' if 'final_score' in frame else 'evidence_score'; top=frame.sort_values(score_col,ascending=False).head(20)[['candidate']]
 b=pd.read_csv(bio,sep='\t').merge(top,on='candidate',how='inner'); b.to_csv(out/'top20_biological_audit.tsv',sep='\t',index=False)
 print(f'Post-ranking biological summary: {len(b)} candidates')
