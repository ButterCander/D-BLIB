import sys, json, pandas as pd, numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
GS={'PI':'PI103','PD':'PD0325901','OXO':'5Z','PIOXO':'PI103+5Z','PIPD':'PI103+PD'}
def load(p): return pd.read_csv(base+p,sep='\t').set_index('task_id')['score']
T={k:load(f'out/tidex_tsir_{k}_24h.tsv') for k in GS}
G={k:load(f'out/tide_p10000_{v}.tsv') for k,v in GS.items()}
M=pd.DataFrame({c:{r:stats.spearmanr(*[pd.concat([G[r],T[c]],axis=1).dropna().iloc[:,i] for i in (0,1)])[0] for r in GS} for c in GS})
print(M.round(2))
d=np.diag(M.values); off=M.values[~np.eye(5,dtype=bool)]
print('matched mean %.2f min %.2f max %.2f; unmatched mean %.2f max %.2f'%(d.mean(),d.min(),d.max(),off.mean(),off.max()))
# TIDE vs naive
tasks=json.load(open(base+'out/task_genes.json'))
for c in ['PI103','PD0325901','5Z','PI103+5Z','PI103+PD']:
    de=pd.read_csv(base+f'out/de_{c}.tsv',sep='\t',index_col=0)
    lf=de['lfc_shrunk'].where(de['padj'].fillna(1)<.05,0.0)
    naive=pd.Series({t:lf.reindex(g).dropna().mean() for t,g in tasks.items() if lf.reindex(g).dropna().size>=1})
    tide=pd.read_csv(base+f'out/tide_p10000_{c}.tsv',sep='\t').set_index('task_id')['score']
    j=pd.concat([tide,naive],axis=1,keys=['tide','naive']).dropna()
    print(c,len(j),'Pearson %.2f'%stats.pearsonr(j.tide,j.naive)[0])
