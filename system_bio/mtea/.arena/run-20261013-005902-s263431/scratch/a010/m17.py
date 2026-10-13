import sys, pickle, pandas as pd, numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
b=base+'.arena/run-20261013-005902-s263431/scratch/a010/'
d=pickle.load(open(b+'m4.pkl','rb'))
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); 
    r=p[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
for k in ['MEKi','TAKi','PI3Ki-TAKi','PI3Ki-MEKi','PI3Ki']:
    v=d[k].copy(); v['q']=bh(np.minimum(1,2*v.p.values))
    print(k, {int(t):(round(float(v[v.id==t].p.iloc[0]),4), round(float(v[v.id==t].q.iloc[0]),2)) for t in (138,140,141)})
# TDO2 check
x=pd.ExcelFile(base+'si/MOESM2.xlsx')
for sh in x.sheet_names:
    df=pd.read_excel(x,sh)
    r=df[df.geneSymbol.isin(['TDO2','IDO1','IDO2'])][['condition','geneID','geneSymbol','log2FoldChange','padj']]
    print(r.to_string(index=False))
hg=pd.read_csv(base+'mteapy/src/mteapy/data/HumanGEM_genes.tsv',sep='\t')
print(hg[hg.geneSymbol.isin(['TDO2','IDO1','IDO2'])].to_string())
# cross-correlation to settle sheet naming
mine={}
for k,f in {'MEKi':'PD0325901','TAKi':'5Z'}.items():
    mine[k]=pd.read_csv(base+f'out/tide_all_p10000_{f}.tsv',sep='\t').set_index('task_id')['score']
for sk in ['MEKi','TAKi']:
    a=d[sk].copy(); a['task_id']='X'+a.id.astype(int).astype(str); a=a.set_index('task_id')['score']
    for mk in ['MEKi','TAKi']:
        j=pd.concat([a,mine[mk]],axis=1).dropna()
        print('sheet',sk,'vs mine',mk,'rho %.2f'%stats.spearmanr(j.iloc[:,0],j.iloc[:,1])[0])
