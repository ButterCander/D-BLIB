import numpy as np, pandas as pd
from scipy import stats
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
C=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0); C=C.loc[:,~C.columns.duplicated()]
ind={'CAPAN2':1.06,'CFPAC1':0.68,'PANC1':0.725,'MIAPACA2':0.155,'KP4':0.18,'ASPC1':-0.10,'SW1990':-0.17}
ks=[k for k in ind if n2i[k] in C.index]; print('lines with CRISPR',ks)
for g in ['SCAP','MBTPS1','SREBF2','MBTPS2','HMGCS1','HMGCR']:
    v=[C.loc[n2i[k],g] for k in ks]; r,p=stats.spearmanr([ind[k] for k in ks],v); print(g,'rho %.2f p %.3f'%(r,p), dict(zip(ks,np.round(v,2))))
# what other genes in all CRISPR correlate with inducibility at |rho|>=0.89 among n=6: how many by chance?
X=C.loc[[n2i[k] for k in ks]].dropna(axis=1)
y=np.array([ind[k] for k in ks])
rhos=np.array([stats.spearmanr(y,X[c].values)[0] for c in X.columns[:6000]])
print('genes tested',len(rhos),'frac |rho|>=0.89: %.3f'%np.mean(np.abs(rhos)>=0.89),' frac <=-0.94 %.3f'%np.mean(rhos<=-0.94))
