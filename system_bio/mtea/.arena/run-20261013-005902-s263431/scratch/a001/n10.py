import numpy as np, pandas as pd
from scipy import stats
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
d=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
d['sym']=gi.reindex(d.index)['gene_name']
d=d.dropna(subset=['sym']); d=d[~d.sym.duplicated()].set_index('sym').drop(columns=[])
d=d.loc[d.sum(axis=1)>=10]
x=d.values.astype(float); pos=(x>0).all(axis=1)
sf=np.median(x[pos]/np.exp(np.log(x[pos]).mean(axis=1,keepdims=True)),axis=0)
Y=pd.DataFrame(np.log2(x/sf+1),index=d.index,columns=d.columns)
def lfc(a,b):
    A=Y[[c for c in Y.columns if c.startswith(a)]]; B=Y[[c for c in Y.columns if c.startswith(b)]]
    lf=B.mean(axis=1)-A.mean(axis=1)
    t,p=stats.ttest_ind(B,A,axis=1,equal_var=False)
    return lf,pd.Series(p,index=Y.index)
for name,(a,b) in {'KP4':('KP4-DMSO','KP4-T100'),'MiaPaCa2':('Mia-DMSO','Mia-T100')}.items():
    lf,p=lfc(a,b)
    g=[s for s in SREBP if s in lf.index]
    print(name,'program n=%d mean %.2f'%(len(g),lf.reindex(g).mean()),'| up',int((lf.reindex(g)>0).sum()),'/',len(g),
          '| HMGCS1 %+.2f (p=%.3f) HMGCR %+.2f MVD %+.2f MVK %+.2f IDI1 %+.2f INSIG1 %+.2f LDLR %+.2f'%tuple([lf['HMGCS1'],p['HMGCS1'],lf['HMGCR'],lf['MVD'],lf['MVK'],lf['IDI1'],lf['INSIG1'],lf['LDLR']]))
    # random 24-gene sets expressed
    rng=np.random.default_rng(5); pool=lf.index.values
    r=np.array([lf.reindex(rng.choice(pool,len(g),replace=False)).mean() for _ in range(5000)])
    print('   random-set percentile of program mean %.3f ; random set SD %.3f'%((r<lf.reindex(g).mean()).mean(), r.std()))
