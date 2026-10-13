import sys, pandas as pd, numpy as np, gzip
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
d=pd.read_csv(base+'pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t')
print(d.shape, d.columns.tolist())
d=d.drop_duplicates('gene_name').set_index('gene_name')
cols=[c for c in d.columns if c not in ('gene_id','gene_name')]
pre=sorted(set(c.rsplit('_',1)[0] for c in cols))
print(pre)
for p in pre: pass
# guess: <LINE>_<M or NT>
lines=sorted(set(p.split('_')[0] for p in pre))
print(lines)
for l in lines:
    treat=[p for p in pre if p.startswith(l+'_') and not p.endswith('_NT')]
    nt=[p for p in pre if p.startswith(l+'_') and p.endswith('_NT')]
    for t in treat:
        lf=np.log2(d[[f'{t}_{i}' for i in (1,2,3)]]+1).mean(axis=1)-np.log2(d[[f'{nt[0]}_{i}' for i in (1,2,3)]]+1).mean(axis=1)
        pr=lf.reindex([g for g in SREBP if g in lf.index]).dropna()
        print(l,t,'vs',nt[0],'program mean %.2f n=%d up=%d'%(pr.mean(),len(pr),(pr>0).sum()), {g:round(lf.get(g,np.nan),2) for g in ['HMGCS1','HMGCR','MVK','MVD']})
