#!/usr/bin/env python3
from pathlib import Path
import numpy as np,pandas as pd
from scipy.stats import spearmanr,pearsonr
from statsmodels.api import OLS,add_constant
from statsmodels.stats.multitest import multipletests
R=Path(__file__).resolve().parents[2]; O=R/'results/phase2_validation'; seed=17; rng=np.random.default_rng(seed)
F=pd.read_csv(O/'top3_mouse_features.tsv',sep='\t'); Q=pd.read_csv(O/'top3_mouse_feature_qc.tsv',sep='\t'); C=['NK → Spp1 → S1pr1 → CD8','T_cell → Cd28 → Cd86 → CD8'];
for c in C:
 q=Q[Q.candidate.eq(c)]; assert len(q)==1 and q.iloc[0].candidate_feature_status=='READY_FOR_STATISTICAL_VALIDATION' and q.iloc[0].n_complete_cases>=6
assert set(Q.loc[Q.candidate.eq('T_cell → Lamc1 → Itga2_Itgb1 → CD4'),'candidate_feature_status'])=={'NOT_ASSESSABLE'}
D=F[F.candidate.isin(C)].copy(); D.to_csv(O/'phase2a_analysis_table.tsv',sep='\t',index=False); rows=[]
for c,g in D.groupby('candidate'):
 x=pd.to_numeric(g.candidate_axis,errors='coerce'); y=pd.to_numeric(g.CD8_activation_consensus,errors='coerce'); ok=x.notna()&y.notna(); x=x[ok].to_numpy(); y=y[ok].to_numpy(); rho,p=spearmanr(x,y); pr,pp=pearsonr(x,y); fit=OLS(y,add_constant(x)).fit(); boots=[]
 for _ in range(5000):
  ix=rng.integers(0,len(x),len(x)); z=spearmanr(x[ix],y[ix]).statistic
  if np.isfinite(z): boots.append(z)
 ci=np.asarray(fit.conf_int())[1]
 rows.append({'candidate':c,'n_mice':len(x),'predictor':'candidate_axis','endpoint':'CD8_activation_consensus','spearman_rho':rho,'spearman_p':p,'spearman_ci_low':np.quantile(boots,.025),'spearman_ci_high':np.quantile(boots,.975),'n_valid_bootstraps':len(boots),'pearson_r':pr,'pearson_p':pp,'ols_beta':fit.params[1],'ols_se':fit.bse[1],'ols_ci_low':ci[0],'ols_ci_high':ci[1],'ols_p':fit.pvalues[1],'ols_r2':fit.rsquared,'direction':'positive' if rho>0 else 'negative' if rho<0 else 'undetermined','association_status':'ASSOCIATION_ONLY'})
A=pd.DataFrame(rows)
A=pd.concat([A,pd.DataFrame([{'candidate':'T_cell → Lamc1 → Itga2_Itgb1 → CD4','n_mice':0,'predictor':'candidate_axis','endpoint':'CD8_activation_consensus','association_status':'NOT_ASSESSABLE'}])],ignore_index=True)
A['spearman_fdr']=np.nan
valid=A.spearman_p.notna()
if valid.any(): A.loc[valid,'spearman_fdr']=multipletests(A.loc[valid,'spearman_p'],method='fdr_bh')[1]
A.to_csv(O/'top3_mouse_level_associations.tsv',sep='\t',index=False)
adj=[]; cond=[]; conf=[]; lomo=[]; summ=[]
for c,g in D.groupby('candidate'):
 x=pd.to_numeric(g.candidate_axis,errors='coerce'); y=pd.to_numeric(g.CD8_activation_consensus,errors='coerce'); ok=x.notna()&y.notna(); x=x[ok].to_numpy(); y=y[ok].to_numpy()
 fit=OLS(y,add_constant(x)).fit(); adj.append({'candidate':c,'model_id':'RAW_OLS','formula':'endpoint ~ candidate_axis','n_mice':len(x),'candidate_axis_beta':fit.params[1],'candidate_axis_se':fit.bse[1],'candidate_axis_ci_low':np.asarray(fit.conf_int())[1,0],'candidate_axis_ci_high':np.asarray(fit.conf_int())[1,1],'candidate_axis_p':fit.pvalues[1],'r2':fit.rsquared,'adjusted_r2':fit.rsquared_adj,'status':'ESTIMABLE','reason':''})
 for i in range(len(x)):
  xx=np.delete(x,i); yy=np.delete(y,i); rr=spearmanr(xx,yy).statistic; ff=OLS(yy,add_constant(xx)).fit(); lomo.append({'candidate':c,'removed_mouse':g.loc[ok].iloc[i].mouse_id,'n_remaining':len(xx),'spearman_rho':rr,'ols_beta':ff.params[1],'direction':'positive' if rr>0 else 'negative' if rr<0 else 'undetermined'})
 z=np.array([r['spearman_rho'] for r in lomo if r['candidate']==c]); b=np.array([r['ols_beta'] for r in lomo if r['candidate']==c]); raw=np.sign(spearmanr(x,y).statistic); frac=np.mean(np.sign(z)==raw); summ.append({'candidate':c,'same_direction_fraction':frac,'rho_min':z.min(),'rho_max':z.max(),'rho_median':np.median(z),'beta_min':b.min(),'beta_max':b.max(),'beta_median':np.median(b),'single_mouse_flip':bool(np.any(np.sign(z)!=raw)),'status':'ROBUST' if frac>=.8 and not np.any(np.sign(z)!=raw) else 'MODERATE' if frac>=.6 else 'FRAGILE'})
 cond.append({'candidate':c,'analysis_type':'OVERALL','condition_or_model':'all','n_mice':len(x),'rho':spearmanr(x,y).statistic,'p':spearmanr(x,y).pvalue,'direction':'positive' if raw>0 else 'negative','status':'DESCRIPTIVE','notes':''})
 conf.append({'candidate':c,'confounder_type':'composition','test_or_model':'not precomputed','status':'NOT_ESTIMABLE','notes':'Covariate unavailable'})
pd.DataFrame(adj).to_csv(O/'top3_adjusted_models.tsv',sep='\t',index=False); pd.DataFrame(cond).to_csv(O/'top3_condition_robustness.tsv',sep='\t',index=False); pd.DataFrame(conf).to_csv(O/'top3_confounding_checks.tsv',sep='\t',index=False); pd.DataFrame(lomo).to_csv(O/'top3_leave_one_mouse_out.tsv',sep='\t',index=False); pd.DataFrame(summ).to_csv(O/'top3_lomo_summary.tsv',sep='\t',index=False)
print(f'Phase 2A statistical validation complete for {len(C)} candidates')
