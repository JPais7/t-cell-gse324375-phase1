#!/usr/bin/env python3
from pathlib import Path
import pandas as pd,anndata as ad
R=Path(__file__).resolve().parents[2];O=R/'results/phase2_validation';O.mkdir(exist_ok=True)
C=['NK → Spp1 → S1pr1 → CD8','T_cell → Cd28 → Cd86 → CD8','T_cell → Lamc1 → Itga2_Itgb1 → CD4']
a=ad.read_h5ad(R/'results/activation/activation_strata_pseudobulk_counts.h5ad');o=a.obs
rows=[]
for c in C:
 for m in sorted(o.mouse_id.unique()): rows.append({'candidate':c,'mouse_id':m,'condition':'NOT_AVAILABLE','time':'NOT_AVAILABLE','checkpoint_blockade':'NOT_AVAILABLE','library_or_batch':'NOT_AVAILABLE'})
pd.DataFrame(rows).to_csv(O/'top3_mouse_features.tsv',sep='\t',index=False)
