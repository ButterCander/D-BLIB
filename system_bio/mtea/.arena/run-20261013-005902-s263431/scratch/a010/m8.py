import sys, pickle, json, pandas as pd, numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
b=base+'.arena/run-20261013-005902-s263431/scratch/a010/'
d=pickle.load(open(b+'m4.pkl','rb'))
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); 
    r=p[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
mp={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
print('A) mine (all-genes null, 10000 perms) vs authors')
for k,f in mp.items():
    m=pd.read_csv(base+f'out/tide_all_p10000_{f}.tsv',sep='\t')
    m2=pd.read_csv(base+f'out/tide_p10000_{f}.tsv',sep='\t')
    a=d[k].copy(); a['task_id']='X'+a.id.astype(int).astype(str)
    j=a.merge(m,on='task_id',suffixes=('_a','_m'))
    r=np.corrcoef(j.score_a,j.score_m)[0,1]; sp=stats.spearmanr(j.score_a,j.score_m)[0]
    q=bh(np.minimum(1,2*m.pvalue.values)); q2=bh(np.minimum(1,2*m2.pvalue.values))
    print(f'{k:11} n={len(j)} r={r:.3f} rho={sp:.3f} | mine all-genes p<.025: {(m.pvalue<0.025).sum()} BH: {(q<0.05).sum()} | model-gene null p<.025 {(m2.pvalue<0.025).sum()} BH {(q2<0.05).sum()}')
print()
print('B) E-Flux grid')
g=json.load(open(base+'out/eflux_grid.json'))
print(type(g), len(g)); 
print(json.dumps(g[0])[:600] if isinstance(g,list) else list(g.keys())[:10])
vals=[v for x in g for v in x['pred'].values()]
print('n values',len(vals),'min',min(vals),'max',max(vals),'n settings',len(g))
