import numpy as np, pandas as pd, itertools
from scipy import stats
P='D:/Cander/School/3-2/system_bio/mtea/'
model=pd.read_csv(P+'external/Model_24Q4.csv').set_index('ModelID')
name2id={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
IND={'CAPAN2':[1.06],'CFPAC1':[1.31,0.05],'PANC1':[1.12,0.33],'MIAPACA2':[0.33,-0.02],'KP4':[0.18],'ASPC1':[-0.10],'SW1990':[-0.17]}
ind=pd.Series({k:np.mean(v) for k,v in IND.items()}).sort_values(ascending=False)
E=pd.read_csv(P+'external/expr_pancreas.tsv',sep='\t',index_col=0)
sel=[name2id[k] for k in ind.index]
print(ind.round(3).to_dict())
trunk=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS']
bt=E.loc[sel,trunk].mean(axis=1)
print('rho trunk',stats.spearmanr(ind.values,bt.values))
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
G=[g for g in SREBP if g in E.columns]
print('program genes present',len(G))
rho_g={g:stats.spearmanr(ind.values,E.loc[sel,g].values)[0] for g in G}
print({g:round(v,2) for g,v in rho_g.items()})
# random 6-gene subsets from the 24 program genes
combs=list(itertools.combinations(G,6))
print('n subsets',len(combs))
rs=np.array([stats.spearmanr(ind.values,E.loc[sel,list(c)].mean(axis=1).values)[0] for c in combs[::max(1,len(combs)//20000)]])
print('subsets n',len(rs),'frac rho>=0.96',(rs>=0.96).mean(),'frac |rho|>=0.96',(np.abs(rs)>=0.96).mean(),'median',np.median(rs),'5-95%',np.percentile(rs,[5,95]))
# random 6-gene sets from expressed genes (their null)
rng=np.random.default_rng(3)
pool=[g for g in E.columns if E[g].median()>2]
rr=np.array([stats.spearmanr(ind.values,E.loc[sel,list(rng.choice(pool,6,replace=False))].mean(axis=1).values)[0] for _ in range(3000)])
print('random expressed gene sets: frac |rho|>=0.96',(np.abs(rr)>=0.96).mean(),'frac rho>=0.96',(rr>=0.96).mean())
# the 4 trametinib-24h-only lines (single study): baseline vs induction
ind4={'CAPAN2':1.06,'CFPAC1':1.31,'MIAPACA2':0.33,'PANC1':1.12}
ids4=[name2id[k] for k in ind4]
print('GSE238198 only (n=4): rho(trunk baseline, 24h program)',stats.spearmanr(list(ind4.values()),E.loc[ids4,trunk].mean(axis=1).values))
# exact permutation p for n=7 rho=0.96
from itertools import permutations
x=np.arange(7); r0=stats.spearmanr(ind.values,bt.values)[0]
cnt=0;tot=0
for perm in permutations(range(7)):
    tot+=1
    if abs(stats.spearmanr(x,perm)[0])>=abs(r0)-1e-9: cnt+=1
print('exact two-sided perm p for n=7 |rho|>=',round(r0,3),cnt/tot)
