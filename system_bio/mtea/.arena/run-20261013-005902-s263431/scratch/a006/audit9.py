import pandas as pd, numpy as np
from scipy import stats
B='D:/Cander/School/3-2/system_bio/mtea/'
model=pd.read_csv(B+'external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv(B+'external/expr_pancreas.tsv',sep='\t',index_col=0)
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
names={'CAPAN2':'CAPAN2','PANC1':'PANC1','CFPAC1':'CFPAC1','KP4':'KP4','MIAPACA2':'MIAPACA2','ASPC1':'ASPC1','SW1990':'SW1990'}
for k in names: print(k, k in n2i, n2i.get(k))
base={k:E.loc[n2i[k],TR].mean() for k in names}
print({k:round(v,2) for k,v in base.items()})
ind={"CAPAN2":1.06,"PANC1":0.725,"CFPAC1":0.68,"KP4":0.18,"MIAPACA2":0.155,"ASPC1":-0.10,"SW1990":-0.17}
x=[base[k] for k in ind]; y=list(ind.values())
print('pooled n=7 rho %.2f p=%.4f'%stats.spearmanr(x,y))
# within GSE238198 trametinib 24h only (n=4)
tr={'CAPAN2':1.06,'CFPAC1':1.31,'MIAPACA2':0.33,'PANC1':1.12}
r,p=stats.spearmanr([base[k] for k in tr],list(tr.values())); print('trametinib-24h-only n=4 rho %.2f p=%.2f'%(r,p))
# trametinib-only lines incl KP4 48h (n=5): 
tr5=dict(tr); tr5['KP4']=0.18
r,p=stats.spearmanr([base[k] for k in tr5],list(tr5.values())); print('+KP4 48h n=5 rho %.2f'%r)
# MAPK-inhibitor class only? (all trametinib incl MiaPaCa-2 48h composite)
# leave out the two KRAS-inhibitor-only lines
sub={k:v for k,v in ind.items() if k not in ('ASPC1','SW1990')}
r,p=stats.spearmanr([base[k] for k in sub],list(sub.values())); print('without AsPC-1,SW1990 (n=5) rho %.2f p=%.2f'%(r,p))
# per-line baseline ranking among all pancreatic lines in expression file
allb=E[TR].mean(axis=1)
print('baseline percentile among %d pancreas models:'%len(allb),{k:round((allb<base[k]).mean(),2) for k in names})
