import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore')
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
C=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0); C=C.loc[:,~C.columns.duplicated()]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
comp={'CAPAN2':1.06,'CFPAC1':(1.31+0.05)/2,'PANC1':(1.12+0.33)/2,'MIAPACA2':(0.33-0.02)/2,'KP4':0.18,'ASPC1':-0.10,'SW1990':-0.17}
ks=list(comp); ids=[n2i[k] for k in ks]; y=np.array([comp[k] for k in ks])
bt=E.loc[ids,TR].mean(axis=1).values
print('obs',stats.spearmanr(bt,y)[0])
rng=np.random.default_rng(7)
pool=[g for g in E.columns if E[g].median()>2]
rs=np.array([stats.spearmanr(E.loc[ids,list(rng.choice(pool,6,replace=False))].mean(axis=1).values,y)[0] for _ in range(5000)])
print('random 6-gene |rho|>=0.964:',(np.abs(rs)>=0.964).mean(),' >=:',(rs>=0.964).mean(), '95th |rho|',np.percentile(np.abs(rs),95))
# random SREBP-program genes? skip
# dependency
for g in ['SCAP','MBTPS1','MBTPS2','SREBF2','HMGCS1']:
    ok=[k for k in ks if n2i[k] in C.index]
    v=C.loc[[n2i[k] for k in ok],g].values; r,p=stats.spearmanr(np.array([comp[k] for k in ok]),v)
    print(g,len(ok),round(r,2),round(p,3))
# baseline of individual genes
for g in ['INSIG1','MBTPS1','SQLE','MVK','HMGCR','MVD','SPRY2','DUSP6','ETV4']:
    r,p=stats.spearmanr(E.loc[ids,g].values,y); print(g,round(r,2),round(p,3))
