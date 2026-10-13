import re,sys
src=open('eflux_grid.py',encoding='utf8').read().split('MEAS = {')[0]
src=src.replace('VMAX = [float(v) for v in sys.argv[1:]] or [10.0]','VMAX=[10.0]')
exec(compile(src,'eflux_grid_head','exec'))
import numpy as np, pandas as pd
ctrl_cols=[c for i,c in enumerate(E.columns) if i%6==0]
col=ctrl_cols[0]
def solve_kd(genes, factor, vmax, floor):
    global RE
    e=E[col].copy()
    ids=[gm[g] for g in genes if g in gm.index]
    for i in ids:
        if i in e.index: e.loc[i]=e.loc[i]*factor
    ed=e.to_dict()
    re=pd.Series({rid: ev(node, ed) for rid,node in gpr.items()})
    RE['__kd']=re
    return solve('__kd', vmax, floor)
base={}
for vmax,floor in [(1.0,0.05),(3.0,0.05),(3.0,0.2),(10.0,0.05)]:
    b=solve(col,vmax,floor); base[(vmax,floor)]=b
    out=[f'vmax={vmax} floor={floor} base={b:.4f}']
    for name,gs in [('RRM1+RRM2',['RRM1','RRM2']),('TYMS',['TYMS']),('DHODH',['DHODH']),('PHGDH',['PHGDH']),('HMGCR',['HMGCR']),('GAPDH',['GAPDH']),('ASNS',['ASNS'])]:
        v=solve_kd(gs,0.1,vmax,floor)
        out.append(f'{name}:{100*(v/b-1):+.1f}%')
    print(' | '.join(out),flush=True)
# random genes 50%
rng=np.random.default_rng(1)
mg=sorted({g for r in model.reactions for g in [x.id for x in r.genes]})
vmax,floor=3.0,0.05
b=base[(vmax,floor)]
v=solve_kd([g for g in []],1,vmax,floor)
# 200 random model genes ×0.5 via ids directly
e=E[col].copy(); ids=[i for i in rng.choice(e.index.values,200,replace=False)]
for i in ids: e.loc[i]=e.loc[i]*0.5
ed=e.to_dict(); RE['__kd']=pd.Series({rid: ev(node, ed) for rid,node in gpr.items()})
print('200 random genes x0.5:',f'{100*(solve("__kd",vmax,floor)/b-1):+.1f}%')
