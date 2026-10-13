import sys, time, numpy as np, pandas as pd
sys.path.insert(0,'mteapy/src')
from cobra.io import read_sbml_model
from mteapy.tide import calculate_TIDE_scores
D='mteapy/src/mteapy/data'
t0=time.time()
hg=pd.read_csv(f'{D}/HumanGEM_genes.tsv',sep='\t').drop_duplicates('geneSymbol')
sym2gem=hg.set_index('geneSymbol')['geneID']; gem_ids=set(hg.geneID)
ts=pd.read_csv(f'{D}/task_structure_matrix.tsv',sep='\t',index_col=0).astype(bool)
model=read_sbml_model(f'{D}/HumanGEM.xml.gz')
gpr={r.id:r.gpr for r in model.reactions if r.id in ts.index}
print('model loaded',time.time()-t0,flush=True)
x=pd.ExcelFile('si/MOESM2.xlsx'); auth_map={}
for sh in x.sheet_names:
    d=x.parse(sh)
    for s,g in zip(d.geneSymbol,d.geneID): auth_map.setdefault(s,set()).add(g)
conflict={s:g for s,g in auth_map.items() if s in sym2gem.index and sym2gem[s] not in g}
print('symbols in authors DEG lists',len(auth_map),'conflicting IDs',len(conflict),list(conflict.items())[:5],flush=True)
print('TDO2',auth_map.get('TDO2'),sym2gem.get('TDO2'))
COND={'PI103':'PI3Ki','PD0325901':'MEKi','5Z':'TAKi','PI103+5Z':'PI3Ki-TAKi','PI103+PD':'PI3Ki-MEKi'}
res={}
for c,sh in COND.items():
    r0=pd.read_csv(f'out/de_{c}.tsv',sep='\t',index_col=0); r0['padj']=r0['padj'].fillna(1)
    for fixed in (False,True):
        ids=[]
        for s in r0.index:
            if fixed and s in auth_map:
                g=[i for i in auth_map[s] if i in gem_ids]; ids.append(g[0] if g else None)
            else: ids.append(sym2gem.get(s))
        r=r0.copy(); r['id']=ids; m=r[r.id.notna()].drop_duplicates('id')
        gd=dict(zip(m.id,np.where(m.padj<.05,m.lfc_shrunk,0.0)))
        sc=calculate_TIDE_scores(gd,ts,gpr,'absmax')
        s=pd.Series(sc,index=[int(t[1:]) for t in ts.columns])
        a=pd.read_excel('si/MOESM4.xlsx',sheet_name=sh,header=None,skiprows=2).iloc[:,[0,2]]; a.columns=['id','s']
        a=a[pd.to_numeric(a.id,errors='coerce').notna()]; a=pd.Series(pd.to_numeric(a.s).values,index=a.id.astype(float).astype(int))
        j=pd.concat([a,s],axis=1,keys=['a','m']).dropna()
        res[(c,fixed)]=(np.corrcoef(j.a,j.m)[0,1], j.corr(method='spearman').iloc[0,1], len(j))
        print(f'{c:10} {sh:11} fixed={fixed} r={res[(c,fixed)][0]:.3f} rho={res[(c,fixed)][1]:.3f} n={len(j)}',flush=True)
print('done',time.time()-t0)
