import sys, pandas as pd, numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
model=pd.read_csv(base+'external/Model_24Q4.csv').set_index('ModelID')
print(model.columns.tolist()[:20])
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
print({k:n2i[k] for k in ['CAPAN2','PANC1','CFPAC1','KP4','MIAPACA2','ASPC1','SW1990'] if k in n2i})
E=pd.read_csv(base+'external/expr_pancreas.tsv',sep='\t',index_col=0)
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
ind={"CAPAN2":1.06,"PANC1":0.725,"CFPAC1":0.68,"KP4":0.18,"MIAPACA2":0.155,"ASPC1":-0.10,"SW1990":-0.17}
x=[E.loc[n2i[k],TR].mean() for k in ind]; y=list(ind.values())
for k,xx,yy in zip(ind,x,y): print(k,round(xx,2),yy)
print('spearman all 7',stats.spearmanr(x,y))
# trametinib-only 24h (GSE238198) values: CAPAN2 1.06, CFPAC1 1.31, PANC1 1.12, MIA 0.33
tr={'CAPAN2':1.06,'CFPAC1':1.31,'PANC1':1.12,'MIAPACA2':0.33}
xx=[E.loc[n2i[k],TR].mean() for k in tr]
print('within GSE238198 24h trametinib n=4 spearman',stats.spearmanr(xx,list(tr.values())))
print('values',dict(zip(tr,[round(v,2) for v in xx])))
# individual dataset-specific: MEKi-only lines: CAPAN2, CFPAC1, PANC1, MIA(24h .33), KP4 (48h .18), MIA 48h -0.02
# drug-class split: lines tested with KRAS inhibitor only
print('KRASi-only lines AsPC-1/SW1990 baseline:', [round(E.loc[n2i[k],TR].mean(),2) for k in ['ASPC1','SW1990']])
print('trametinib-24h lines baseline:', {k:round(E.loc[n2i[k],TR].mean(),2) for k in tr})
# pancreas distribution of baseline trunk expression among all 55 lines
allm=E[TR].mean(axis=1)
print('55 lines baseline trunk mean: min',allm.min().round(2),'median',allm.median().round(2),'max',allm.max().round(2), 'n',len(allm))
rank=allm.rank(pct=True)
for k in ind: print(k, 'percentile among pancreatic lines', round(rank[n2i[k]],2))
# random 6-gene null as in report
rng=np.random.default_rng(0)
genes=[g for g in E.columns if E[g].notna().all()]
cnt=0; N=3000
for _ in range(N):
    gs=rng.choice(genes,6,replace=False)
    xx=[E.loc[n2i[k],gs].mean() for k in ind]
    r=stats.spearmanr(xx,y)[0]
    if abs(r)>=0.96: cnt+=1
print('random 6-gene sets |rho|>=0.96:',cnt/N)
# among genes within mevalonate/SREBP list: would any other trunk-like set do as well? 
