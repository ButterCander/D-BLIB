import numpy as np, pandas as pd
from scipy import stats
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
print('expr shape',E.shape)
TR=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS']
lines=['CAPAN2','CFPAC1','PANC1','MIAPACA2','KP4','ASPC1','SW1990']
base={l:E.loc[n2i[l],TR].mean() for l in lines}
print({k:round(v,2) for k,v in base.items()})
# inducibility variants (values from the report; recomputed program means from datasets below)
ind_mix={'CAPAN2':1.06,'CFPAC1':(1.31+0.05)/2,'PANC1':(1.12+0.33)/2,'MIAPACA2':(0.33-0.02)/2,'KP4':0.18,'ASPC1':-0.10,'SW1990':-0.17}
x=[base[l] for l in lines]; y=[ind_mix[l] for l in lines]
print('mixed-drug composite rho %.3f p %.4f'%stats.spearmanr(x,y))
# trametinib-only (24h GSE238198 for 4 lines + KP4 48h)
tram={'CAPAN2':1.06,'CFPAC1':1.31,'PANC1':1.12,'MIAPACA2':0.33,'KP4':0.18}
l2=list(tram); print('trametinib-only (n=5) rho %.3f'%stats.spearmanr([base[l] for l in l2],[tram[l] for l in l2])[0])
tram4={k:tram[k] for k in ['CAPAN2','CFPAC1','PANC1','MIAPACA2']}
print('trametinib 24h only (n=4) rho %.3f'%stats.spearmanr([base[l] for l in tram4],[tram4[l] for l in tram4])[0])
# random 6-gene control
rng=np.random.default_rng(1)
pool=[g for g in E.columns if E[g].median()>2]
sel=[n2i[l] for l in lines]
obs=stats.spearmanr(y,E.loc[sel,TR].mean(axis=1).values)[0]
rs=np.array([stats.spearmanr(y,E.loc[sel,list(rng.choice(pool,6,replace=False))].mean(axis=1).values)[0] for _ in range(3000)])
print('obs %.3f; frac random >= obs %.4f; frac |rho|>=|obs| %.4f; pool %d'%(obs,(rs>=obs).mean(),(np.abs(rs)>=abs(obs)).mean(),len(pool)))
# exact permutation p for n=7, rho=.96
ys=np.array(y); xs=np.array(x)
from itertools import permutations
cnt=0;tot=0
for p in permutations(range(7)):
    r=stats.spearmanr(xs,ys[list(p)])[0]; tot+=1; cnt+= (abs(r)>=abs(obs)-1e-12)
print('exact permutation p (two-sided) %.5f of %d'%(cnt/tot,tot))
# baseline levels in normalized counts for pseudocount concern
for l in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    d=pd.read_csv(f'pdac/GSE238198_Pdac.{l}.genes.normalizedCounts.txt.gz',sep='\t',index_col=0)
    nt=d[[c for c in d.columns if '_NT_' in c]].mean(axis=1)
    print(l,'NT normalized counts trunk:',nt.reindex(TR).round(0).to_dict())
