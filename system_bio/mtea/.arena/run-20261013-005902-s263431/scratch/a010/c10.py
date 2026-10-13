import pandas as pd, numpy as np, itertools, warnings
from scipy import stats
warnings.filterwarnings('ignore')
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
lines=['CAPAN2','CFPAC1','PANC1','MIAPACA2','KP4','ASPC1','SW1990']
for l in lines: print(l, l in n2i, n2i.get(l))
bt={l:E.loc[n2i[l],TR].mean() for l in lines}
print({k:round(v,2) for k,v in bt.items()})
def rho(ind):
    ks=list(ind); return stats.spearmanr([bt[k] for k in ks],[ind[k] for k in ks])
comp={'CAPAN2':1.06,'CFPAC1':(1.31+0.05)/2,'PANC1':(1.12+0.33)/2,'MIAPACA2':(0.33-0.02)/2,'KP4':0.18,'ASPC1':-0.10,'SW1990':-0.17}
print('composite n=7',rho(comp))
tram24={'CAPAN2':1.06,'CFPAC1':1.31,'PANC1':1.12,'MIAPACA2':0.33}
print('trametinib 24h only n=4',rho(tram24))
tram_all={'CAPAN2':1.06,'CFPAC1':1.31,'PANC1':1.12,'MIAPACA2':(0.33-0.02)/2,'KP4':0.18}
print('trametinib 24h + 48h n=5',rho(tram_all))
mrtx={'ASPC1':-0.10,'PANC1':0.33,'CFPAC1':0.05,'SW1990':-0.17}
print('KRASi only n=4',rho(mrtx))
# drug-class-adjusted: subtract the mean of the drug type? lines with both drugs: PANC1 (1.12 vs .33), CFPAC1 (1.31 vs .05) -> KRASi gives lower by ~0.9-1.3
d_shift=np.mean([1.12-0.33,1.31-0.05])
adj=dict(comp)
adj['ASPC1']=-0.10+d_shift; adj['SW1990']=-0.17+d_shift
print('shift of KRASi-only lines by mean tram-minus-KRASi difference',round(d_shift,2),rho(adj))
# if AsPC-1 and SW1990 had trametinib induction equal to the weakest tram line (MiaPaCa2 24h 0.33)
alt=dict(comp); alt['ASPC1']=0.33; alt['SW1990']=0.33
print('if both ASPC1,SW1990 = 0.33 under tram',rho(alt))
# only lines that have trametinib data, using only trametinib datasets, 24h: n=4; baseline trunk
# sensitivity to drop of single genes
for g in TR:
    t=[x for x in TR if x!=g]
    b2={l:E.loc[n2i[l],t].mean() for l in lines}
    print('drop',g, round(stats.spearmanr([b2[k] for k in comp],[comp[k] for k in comp])[0],2))
# exact permutation p for n=7 rho=0.964
x=np.array([bt[k] for k in comp]); y=np.array([comp[k] for k in comp])
r0=stats.spearmanr(x,y)[0]
cnt=0;tot=0
for p in itertools.permutations(range(7)):
    tot+=1
    if stats.spearmanr(x,y[list(p)])[0]>=r0-1e-9: cnt+=1
print('exact one-sided perm p',cnt/tot,'rho',r0)
