import numpy as np, pandas as pd, json
from scipy import stats
B='D:/Cander/School/3-2/system_bio/mtea/'
def load(p): return pd.read_csv(B+p,sep='\t').set_index('task_id')['score']
GS={"PI":"PI103","PD":"PD0325901","OXO":"5Z","PIOXO":"PI103+5Z","PIPD":"PI103+PD"}
T={k:load('out/tidex_tsir_%s_24h.tsv'%k) for k in GS}
G={k:load('out/tide_p10000_%s.tsv'%v) for k,v in GS.items()}
M=pd.DataFrame({c:{r:stats.spearmanr(*pd.concat([G[r],T[c]],axis=1).dropna().values.T)[0] for r in GS} for c in GS})
print(M.round(2))
d=np.diag(M.values); off=M.values[~np.eye(5,dtype=bool)]
print('matched mean %.2f min %.2f max %.2f; unmatched mean %.2f max %.2f'%(d.mean(),d.min(),d.max(),off.mean(),off.max()))
# naive vs TIDE
tasks=json.load(open(B+'out/task_genes.json'))
for c in ["PI103","PD0325901","5Z","PI103+5Z","PI103+PD"]:
    de=pd.read_csv(B+'out/de_%s.tsv'%c,sep='\t',index_col=0)
    lf=de['lfc_shrunk'].where(de['padj'].fillna(1)<.05,0.0)
    naive=pd.Series({t:lf.reindex(g).dropna().mean() for t,g in tasks.items() if lf.reindex(g).dropna().size>=1})
    tide=load('out/tide_p10000_%s.tsv'%c)
    j=pd.concat([tide,naive],axis=1,keys=['t','n']).dropna()
    print(c,len(j),'pearson %.2f'%stats.pearsonr(j.t,j.n)[0])
