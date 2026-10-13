import pandas as pd, numpy as np, sys
from scipy.stats import spearmanr, pearsonr
conds=['5Z','PD0325901','PI103','PI103+5Z','PI103+PD']
T={c:pd.read_csv(f'out/tide_{c}.tsv',sep='\t').set_index('task_id') for c in conds}
E={c:pd.read_csv(f'out/tidee_{c}.tsv',sep='\t').set_index('task_id') for c in conds}
rows=[]
for c in conds:
    j=T[c].join(E[c],lsuffix='_t',rsuffix='_e',how='inner')
    j['cond']=c; rows.append(j)
J=pd.concat(rows)
print('n pairs both defined',J.score_e.notna().sum(),'of',len(J))
Jd=J.dropna(subset=['score_e'])
print('Spearman all pairs (both defined):',round(spearmanr(Jd.score_t,Jd.score_e)[0],3))
sig_t=(J.pvalue_t<0.025); sig_e=(J.pvalue_e<0.025)
print('TIDE sig pairs',sig_t.sum(),'TIDE-E sig pairs',sig_e.sum(),'both',(sig_t&sig_e).sum())
b=J[sig_t&sig_e]
print('Spearman both significant:',round(spearmanr(b.score_t,b.score_e)[0],3),'n',len(b))
# nonzero scores only
nz=Jd[(Jd.score_t!=0)&(Jd.score_e!=0)]
print('nonzero both n',len(nz),'rho',round(spearmanr(nz.score_t,nz.score_e)[0],3))
for c in conds:
    print(c,'sig TIDE',(T[c].pvalue<0.025).sum(),'sig TIDEe',(E[c].pvalue<0.025).sum())
