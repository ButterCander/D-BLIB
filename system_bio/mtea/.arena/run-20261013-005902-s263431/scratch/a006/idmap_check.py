import sys, numpy as np, pandas as pd
sys.path.insert(0,'D:/Cander/School/3-2/system_bio/mtea/mteapy/src')
from scipy import stats
from cobra.io import read_sbml_model
from mteapy.tide import calculate_TIDE_scores
B='D:/Cander/School/3-2/system_bio/mtea/'
D=B+'mteapy/src/mteapy/data'
hg=pd.read_csv(D+'/HumanGEM_genes.tsv',sep='\t').drop_duplicates('geneSymbol')
sym2gem=hg.set_index('geneSymbol')['geneID']; gem_ids=set(hg.geneID)
ts=pd.read_csv(D+'/task_structure_matrix.tsv',sep='\t',index_col=0).astype(bool)
model=read_sbml_model(D+'/HumanGEM.xml.gz')
gpr={r.id:r.gpr for r in model.reactions if r.id in ts.index}
x=pd.ExcelFile(B+'si/MOESM2.xlsx'); auth={}
for sh in x.sheet_names:
    d=x.parse(sh)
    for s,g in zip(d.geneSymbol,d.geneID): auth.setdefault(s,set()).add(g)
conf={s:g for s,g in auth.items() if s in sym2gem.index and sym2gem[s] not in g}
print('symbols in authors DEG lists:',len(auth),' with ID != Human-GEM ID:',len(conf))
print('TDO2 authors ID:',auth.get('TDO2'),' HumanGEM ID:',sym2gem.get('TDO2'))
COND={'PI103':'PI3Ki','PD0325901':'MEKi','5Z':'TAKi','PI103+5Z':'PI3Ki-TAKi','PI103+PD':'PI3Ki-MEKi'}
def build(c,fixed):
    r0=pd.read_csv(B+'out/de_%s.tsv'%c,sep='\t',index_col=0); r0['padj']=r0['padj'].fillna(1)
    ids=[]
    for s in r0.index:
        if fixed and s in auth:
            g=[i for i in auth[s] if i in gem_ids]; ids.append(g[0] if g else None)
        else: ids.append(sym2gem.get(s))
    r0['id']=ids; m=r0[r0.id.notna()].drop_duplicates('id')
    return dict(zip(m.id.values,np.where(m.padj<.05,m.lfc_shrunk,0.0)))
for c,sh in COND.items():
    a=pd.read_excel(B+'si/MOESM4.xlsx',sheet_name=sh,header=None).iloc[2:,[0,2]]; a.columns=['id','s']; a['s']=pd.to_numeric(a.s); a['tid']='X'+a.id.astype(int).astype(str)
    out=[]
    for fixed in (False,True):
        sc=calculate_TIDE_scores(build(c,fixed),ts,gpr,'absmax')
        s=pd.Series(sc,index=ts.columns)
        j=a.assign(m=s.reindex(a.tid).values)
        out.append(np.corrcoef(j.s,j.m)[0,1])
    print('%-11s r symbol-map %.2f -> ID-fix %.2f'%(sh,out[0],out[1]))
