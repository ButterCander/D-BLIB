import pandas as pd, numpy as np
from scipy import stats
pd.set_option('display.width',250)
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
print(E.shape, [c for c in E.columns[:5]])
TR=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS']
ind={'CAPAN2':1.06,'PANC1':0.725,'CFPAC1':0.68,'KP4':0.18,'MIAPACA2':0.155,'ASPC1':-0.10,'SW1990':-0.17}
for k in ind: print(k, k in n2i)
x=np.array([E.loc[n2i[k],TR].mean() for k in ind]); y=np.array(list(ind.values()))
print('rho all7',stats.spearmanr(x,y))
print(pd.Series(dict(zip(ind,x.round(2)))))
# within study
tram24={'CAPAN2':1.06,'CFPAC1':1.31,'MIAPACA2':0.33,'PANC1':1.12}
xs=[E.loc[n2i[k],TR].mean() for k in tram24]; print('rho within trametinib 24h (n=4)',stats.spearmanr(xs,list(tram24.values())))
mrtx={'PANC1':0.33,'ASPC1':-0.10,'SW1990':-0.17,'CFPAC1':0.05}
xs=[E.loc[n2i[k],TR].mean() for k in mrtx]; print('rho within MRTX/BI (n=4)',stats.spearmanr(xs,list(mrtx.values())))
# random 6-gene sets
rng=np.random.default_rng(0)
cols=[c for c in E.columns]
expr=E.loc[[n2i[k] for k in ind]]
cnt=0;N=3000
rs=[]
for _ in range(N):
    g=rng.choice(cols,6,replace=False)
    r=stats.spearmanr(expr[g].mean(axis=1).values,y)[0]
    rs.append(r)
rs=np.array(rs)
print('random sets |rho|>=0.96:',np.mean(np.abs(rs)>=0.96),'; rho>=0.96 one-sided',np.mean(rs>=0.96),'; 95th pct |rho|',np.percentile(np.abs(rs),95))
# only expressed (mean>some) genes random sets, to avoid junk genes
hi=[c for c in cols if E[c].mean()>np.median(E.mean())]
rs2=np.array([stats.spearmanr(expr[rng.choice(hi,6,replace=False)].mean(axis=1).values,y)[0] for _ in range(N)])
print('random among above-median-expressed genes: |rho|>=.96',np.mean(np.abs(rs2)>=0.96),'rho>=.96',np.mean(rs2>=0.96))
# leave one out
for k in ind:
    kk=[a for a in ind if a!=k]; print('drop',k, round(stats.spearmanr([E.loc[n2i[a],TR].mean() for a in kk],[ind[a] for a in kk])[0],2))
# the 6 trunk genes themselves are part of the 24-gene program used to define inducibility: circularity check: inducibility excluding trunk genes?
