import sys,time, numpy as np, pandas as pd
sys.path.insert(0,'mteapy/src')
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D='mteapy/src/mteapy/data'
gm=pd.read_csv(f'{D}/HumanGEM_genes.tsv',sep='\t').drop_duplicates('geneSymbol').set_index('geneSymbol')['geneID']
ts=pd.read_csv(f'{D}/task_structure_matrix.tsv',sep='\t',index_col=0)
model=read_sbml_model(f'{D}/HumanGEM.xml.gz')
r0=pd.read_csv('out/de_PI103.tsv',sep='\t',index_col=0); r0['padj']=r0['padj'].fillna(1)
mg=r0.loc[r0.index.intersection(gm.index)].copy(); mg.index=[gm[i] for i in mg.index]
mg['lfc']=np.where(mg.padj<.05,mg.lfc_shrunk,0.0)
print('model genes',len(mg),'nonzero',(mg.lfc!=0).sum(),'tested total',len(r0),flush=True)
au=pd.read_excel('si/MOESM4.xlsx',sheet_name='PI3Ki',header=None,skiprows=2).iloc[:,[0,2,3]]; au.columns=['tid','sa','pa']
au=au[pd.to_numeric(au.tid,errors='coerce').notna()].copy(); au['tid']='X'+au.tid.astype(float).astype(int).astype(str); au['pa']=pd.to_numeric(au.pa); au['sa']=pd.to_numeric(au.sa)
for N in [len(mg),3000,5000,len(r0)]:
    e=pd.DataFrame({'lfc':mg.lfc.values},index=mg.index)
    extra=max(0,N-len(e))
    if extra: e=pd.concat([e,pd.DataFrame({'lfc':np.zeros(extra)},index=[f'dummy{i}' for i in range(extra)])])
    t0=time.time()
    res=compute_TIDE(e,'lfc',ts,model,'absmax',2000,1,True,1)
    nulls={r.task_id:np.array(r.random_score_array.split(';'),dtype=float) for r in res.itertuples()}
    p=[]
    for r in au.itertuples():
        n=nulls[r.tid]; p.append(min((n<=r.sa).mean(),(n>=r.sa).mean()))
    au['pu']=p; q=bh(np.clip(2*au.pu.values,1e-12,1))
    print(f'N={N} nonzero frac {(e.lfc!=0).mean():.3f} authors-scores p<.025: {(au.pu<.025).sum()} BH: {(q<.05).sum()} (authors {(au.pa<.025).sum()} / 59) {time.time()-t0:.0f}s',flush=True)
