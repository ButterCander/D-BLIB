import pandas as pd, numpy as np, json, openpyxl
from scipy import stats
print("=== F. Flobak growth synergy table ===")
f=pd.read_csv('out/flobak_synergy.tsv',sep='\t',index_col=0)
print(f.loc[[18,21],['bliss_excess','bliss_excess_dir','CI_cal','CI_dir','hsa_margin_full','fAB']].round(3))
syn=f[f.call=='syn']; non=f[f.call=='non']
print('n syn',len(syn),'n non',len(non))
print('HSA syn range',syn.hsa_margin_full.min().round(2),syn.hsa_margin_full.max().round(2),'non',non.hsa_margin_full.min().round(2),non.hsa_margin_full.max().round(2))
print('Bliss syn',syn.bliss_excess.min().round(2),syn.bliss_excess.max().round(2),'non',non.bliss_excess.min().round(2),non.bliss_excess.max().round(2))
print('Bliss non mean/sd',non.bliss_excess.mean().round(3),non.bliss_excess.std().round(3))
print('Bliss dir syn',syn.bliss_excess_dir.min().round(2),'non max',non.bliss_excess_dir.max().round(2))
ci=pd.to_numeric(non.CI_cal,errors='coerce'); print('non CI defined',ci.notna().sum(),'of',len(non),'range',ci.min().round(2),ci.max().round(2))
print('syn CI',syn.CI_cal.round(2).tolist())
def auc(a,b): return np.mean([ (x>y)+0.5*(x==y) for x in a for y in b])
print('AUC HSA',auc(syn.hsa_margin_full,non.hsa_margin_full),'Bliss',auc(syn.bliss_excess,non.bliss_excess),'Bliss dir',auc(syn.bliss_excess_dir,non.bliss_excess_dir))
z=(f.loc[18,'bliss_excess']-non.bliss_excess.mean())/non.bliss_excess.std(); z2=(f.loc[21,'bliss_excess']-non.bliss_excess.mean())/non.bliss_excess.std(); print('z',z.round(2),z2.round(2))
# single-drug reproducibility
for drug in ['PI103','PD0325901','5Z']:
    v=[]
    for i,r in f.iterrows():
        if r.A==drug: v.append(r.fA)
        if r.B==drug: v.append(r.fB)
    print(drug,'n',len(v),'min',round(min(v),2),'max',round(max(v),2),'sd',round(np.std(v,ddof=1),3))
# half-dose Bliss threshold test: would the two reference combos pass H2 threshold +0.3?
print('PI3Ki+MEKi excess',round(f.loc[18,'bliss_excess'],2),'PI3Ki+TAKi',round(f.loc[21,'bliss_excess'],2),'threshold 0.3; syn >=0.3:',(syn.bliss_excess>=0.3).sum(),'of',len(syn))
print("-2SD threshold: ", round(non.bliss_excess.mean()+2*non.bliss_excess.std(),2))
print()
print("=== G. E-Flux grid ===")
g=json.load(open('out/eflux_grid.json'))
print(type(g),len(g)); 
vals=[v for x in g for v in x['pred'].values()]
print('min',round(min(vals),1),'max',round(max(vals),1),'n settings',len(g)); print(g[0].keys()); print(g[0])
print()
print("=== H. statin ranks (tidex_statin_*) ===")
for line in ['PANC1','MiaPaca2']:
    for st in ['Fluvastatin','Lovastatin','Simvastatin']:
        for dose in ['High','Low']:
            d=pd.read_csv(f'out/tidex_statin_{line}_{st}_{dose}.tsv',sep='\t').set_index('task_id')
            r=d['score'].dropna().rank(ascending=False)
            ra=d['score'].abs().dropna().rank(ascending=False)
            print(line,st,dose,'ranks X138/X140/X141:',[int(r.get(t,-1)) for t in ('X138','X140','X141')],'n',len(r),'p',[d.loc[t,'pvalue'] if 'pvalue' in d else None for t in ('X138','X140','X141')])
