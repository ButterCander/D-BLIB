import pandas as pd, numpy as np, json
g=json.load(open('out/eflux_grid.json'))
vals=[]
for x in g:
    p=x['pred']; vals+=list(p.values())
    print(x['vmax'],x['floor'],{k:round(v,1) for k,v in p.items()},'rho',round(x['rho'],2),'conc',x.get('conc'))
print(len(g),'settings; min',min(vals),'max',max(vals))
f=pd.read_csv('out/flobak_synergy.tsv',sep='\t',index_col=0)
for d in ['PI103','PD0325901','5Z']:
    v=pd.concat([f[f.A==d].fA,f[f.B==d].fB]); print(d,len(v),v.round(3).tolist(),'min',v.min().round(2),'max',v.max().round(2),'sd',v.std().round(3))
non=f[f.call=='non']
print('non bliss_excess mean sd',non.bliss_excess.mean().round(3),non.bliss_excess.std().round(3),'n',len(non))
print('syn bliss',f[f.call=='syn'].bliss_excess.round(3).tolist(), 'hsa',f[f.call=='syn'].hsa_margin_full.round(3).tolist(),'CI',f[f.call=='syn'].CI_cal.round(3).tolist())
print('non CI',non.CI_cal.dropna().round(2).tolist())
# informative negatives: both single drugs fractional inhib between 0.2 and 1.0
inf=non[(non.fA.between(0.2,1.0))&(non.fB.between(0.2,1.0))]
print(inf[['A','B','fA','fB','fAB','bliss_excess','bliss_excess_dir','hsa_margin_full','CI_cal','CI_dir']].round(3))
trivial=non[(non.A.isin(['CT99021','BIRB0796']))|(non.B.isin(['CT99021','BIRB0796']))]
print('pairs with an inactive partner',len(trivial))
pkf=non[((non.A=='PKF118-310')|(non.B=='PKF118-310'))&~non.index.isin(trivial.index)]
print('pairs with PKF (no inactive)',len(pkf), pkf.index.tolist())
