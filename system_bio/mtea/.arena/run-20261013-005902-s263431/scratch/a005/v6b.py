import pandas as pd, numpy as np
from scipy import stats
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
print(E.shape)
TR=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS']
ind={'CAPAN2':1.06,'PANC1':0.725,'CFPAC1':0.68,'KP4':0.18,'MIAPACA2':0.155,'ASPC1':-0.10,'SW1990':-0.17}
ids=[n2i[k] for k in ind]
M=E.loc[ids].values.astype(float); cols=np.array(E.columns); y=np.array(list(ind.values()))
x=E.loc[ids,TR].mean(axis=1).values
print('rho all7',stats.spearmanr(x,y))
print(dict(zip(ind,x.round(2))))
tram={'CAPAN2':1.06,'CFPAC1':1.31,'MIAPACA2':0.33,'PANC1':1.12}
print('within trametinib24h n=4',stats.spearmanr([E.loc[n2i[k],TR].mean() for k in tram],list(tram.values())))
mr={'PANC1':0.33,'ASPC1':-0.10,'SW1990':-0.17,'CFPAC1':0.05}
print('within KRASi 24h n=4',stats.spearmanr([E.loc[n2i[k],TR].mean() for k in mr],list(mr.values())))
rng=np.random.default_rng(0)
ry=stats.rankdata(y)
def rho(v): 
    r=stats.rankdata(v); return np.corrcoef(r,ry)[0,1]
res=[]
for _ in range(3000):
    j=rng.choice(M.shape[1],6,replace=False); res.append(rho(M[:,j].mean(axis=1)))
res=np.array(res); print('random6 |rho|>=0.96: %.4f; rho>=0.96: %.4f; 95pct |rho| %.2f'%(np.mean(np.abs(res)>=0.96),np.mean(res>=0.96),np.percentile(np.abs(res),95)))
# restrict to genes with mean expr above median (expressed)
mm=M.mean(axis=0); hi=np.where(mm>np.median(mm))[0]
res2=np.array([rho(M[:,rng.choice(hi,6,replace=False)].mean(axis=1)) for _ in range(3000)])
print('random6 among expressed: |rho|>=.96 %.4f rho>=.96 %.4f'%(np.mean(np.abs(res2)>=0.96),np.mean(res2>=0.96)))
for k in ind:
    keep=[i for i,a in enumerate(ind) if a!=k]
    print('drop',k,round(stats.spearmanr(x[keep],y[keep])[0],2))
# SPRY2 / SCAP-type panel and other SREBP genes baseline
for g in ['SPRY2','DUSP6','ETV4','INSIG1','SREBF2','SCAP','LDLR','SQLE','FASN','ACACA']:
    if g in E.columns: print(g, round(stats.spearmanr(E.loc[ids,g].values,y)[0],2))
