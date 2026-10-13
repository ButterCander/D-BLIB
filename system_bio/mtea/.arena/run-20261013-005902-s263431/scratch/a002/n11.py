import numpy as np, pandas as pd
from scipy import stats
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
name2id={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
IND={"CAPAN2":[1.06],"CFPAC1":[1.31,0.05],"PANC1":[1.12,0.33],"MIAPACA2":[0.33,-0.02],"KP4":[0.18],"ASPC1":[-0.10],"SW1990":[-0.17]}
ind=pd.Series({k:np.mean(v) for k,v in IND.items()})
bt=pd.Series({k:E.loc[name2id[k],TR].mean() for k in IND})
print(pd.DataFrame({'induction':ind.round(2),'baseline_trunk_log2TPM':bt.round(2),'n_datasets':{k:len(v) for k,v in IND.items()}}).sort_values('induction',ascending=False))
print('pooled n=7 rho', stats.spearmanr(ind[bt.index],bt)[0].round(3), 'p',stats.spearmanr(ind[bt.index],bt)[1])
# trametinib-24h-only, n=4 (CAPAN2, CFPAC1, MIAPACA2, PANC1)
tram={"CAPAN2":1.06,"CFPAC1":1.31,"MIAPACA2":0.33,"PANC1":1.12}
x=[bt[k] for k in tram]; y=list(tram.values())
print('trametinib 24h only n=4 rho (DepMap baseline)',stats.spearmanr(x,y)[0].round(2), 'pearson',np.corrcoef(x,y)[0,1].round(2))
# MRTX/BI 24h only n=4: AsPC1,PANC1,SW1990,CFPAC1
kr={"ASPC1":-0.10,"PANC1":0.33,"SW1990":-0.17,"CFPAC1":0.05}
x=[bt[k] for k in kr]; y=list(kr.values())
print('KRASi 24h only n=4 rho',stats.spearmanr(x,y)[0].round(2))
# using only first (single) dataset: trametinib wherever available else KRASi
# leave-study: composite by drug
# What if we use trametinib-only inducibility for the 5 lines with trametinib (CAPAN2,CFPAC1,PANC1,MIAPACA2 24h; KP4 48h)
tr5={"CAPAN2":1.06,"CFPAC1":1.31,"MIAPACA2":0.33,"PANC1":1.12,"KP4":0.18}
x=[bt[k] for k in tr5]; y=list(tr5.values()); print('trametinib any time n=5 rho',stats.spearmanr(x,y)[0].round(2))
# exact permutation p for rho=0.96 n=7
import itertools
r=stats.spearmanr(ind[bt.index],bt)[0]
cnt=0;tot=0
rk=stats.rankdata(ind[bt.index].values)
for p in itertools.permutations(range(7)):
    tot+=1
    if abs(stats.spearmanr(rk,np.array(p))[0])>=abs(r)-1e-9: cnt+=1
print('exact two-sided perm p',cnt/tot)
