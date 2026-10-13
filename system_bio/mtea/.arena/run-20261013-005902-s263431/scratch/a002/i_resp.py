import numpy as np, pandas as pd
from scipy import stats
model=pd.read_csv("external/Model_24Q4.csv").set_index("ModelID")
name2id={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=="Pancreas"}
IND={"CAPAN2":[1.06],"CFPAC1":[1.31,0.05],"PANC1":[1.12,0.33],"MIAPACA2":[0.33,-0.02],"KP4":[0.18],"ASPC1":[-0.10],"SW1990":[-0.17]}
ind=pd.Series({k:np.mean(v) for k,v in IND.items()}).sort_values(ascending=False)
E=pd.read_csv("external/expr_pancreas.tsv",sep="\t",index_col=0)
ids={k:name2id[k] for k in IND}
sel=[ids[k] for k in ind.index]
trunk=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
bt=E.loc[sel,trunk].mean(axis=1)
r=stats.spearmanr(ind.values,bt.values)
print("n=",len(ind),"rho",round(r[0],3),"p",round(r[1],4))
print(dict(zip(ind.index,bt.round(2).values)))
# sensitivity: use trametinib-only datasets where available (GSE238198 24h) and drop KRASi-only lines
trm={"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":0.33}
s=pd.Series(trm); b=E.loc[[ids[k] for k in s.index],trunk].mean(axis=1)
print("trametinib-24h only (n=4) rho",round(stats.spearmanr(s.values,b.values)[0],2))
# alt: use only the first-listed dataset value for each line
IND2={"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":0.33,"KP4":0.18,"ASPC1":-0.10,"SW1990":-0.17}
s2=pd.Series(IND2); b2=E.loc[[ids[k] for k in s2.index],trunk].mean(axis=1)
print("first-dataset values (n=7) rho",round(stats.spearmanr(s2.values,b2.values)[0],2))
# exact permutation p for n=7
import itertools
x=stats.rankdata(ind.values); y=stats.rankdata(bt.values)
cnt=0; tot=0; obs=stats.spearmanr(x,y)[0]
for perm in itertools.permutations(range(7)):
    tot+=1
    if abs(stats.spearmanr(x,y[list(perm)])[0])>=abs(obs)-1e-12: cnt+=1
print("exact two-sided perm p",cnt/tot)
# random 6-gene sets
rng=np.random.default_rng(3)
pool=[g for g in E.columns if E[g].median()>2]
rs=np.array([stats.spearmanr(ind.values,E.loc[sel,list(rng.choice(pool,6,replace=False))].mean(axis=1).values)[0] for _ in range(3000)])
print("random 6-gene: frac |rho|>=obs",(np.abs(rs)>=abs(r[0])).mean(),"95th pct |rho|",np.percentile(np.abs(rs),95).round(2))
# SCAP dependency
C=pd.read_csv("external/crispr_pancreas_rows.tsv",sep="\t",index_col=0); C=C.loc[:,~C.columns.duplicated()]
for g in ["SCAP","MBTPS1","MBTPS2","SREBF2","HMGCS1"]:
    ok=[k for k in ind.index if ids[k] in C.index]
    v=C.loc[[ids[k] for k in ok],g].values
    rr=stats.spearmanr(ind[ok].values,v)
    print(g,"n",len(ok),"rho",round(rr[0],2),"p",round(rr[1],3))
print("lines without CRISPR:",[k for k in ind.index if ids[k] not in C.index])
