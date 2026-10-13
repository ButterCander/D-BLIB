import numpy as np, pandas as pd
from scipy import stats
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
print('GSE279749 genes',d.shape[0])
for k,nm in [('A','AsPC-1'),('P','PANC-1'),('C','CFPAC-1'),('S9','SW1990')]:
    nt=np.log2(d[[f'{k}_NT_{i}' for i in (1,2,3)]]+1).mean(axis=1); m=np.log2(d[[f'{k}_M_{i}' for i in (1,2,3)]]+1).mean(axis=1)
    lf=m-nt; g=[s for s in SREBP if s in lf.index]
    hs=[s for s in ['HMGCS1','HMGCR','MVD','IDI1','MVK'] if s in lf.index]
    print(nm,'program n=%d mean %+.2f up %d'%(len(g),lf.reindex(g).mean(),(lf.reindex(g)>0).sum()),{s:round(lf[s],2) for s in hs})
for ln in ['PSN1','SW1990']:
    f=f'pdac/GSE315130_{ln}_cell_line_gene_sample_FPKM_with_symbol.txt.gz'
    d=pd.read_csv(f,sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    cols=d.columns.tolist(); print(ln,[c for c in cols if c.startswith(ln)])
    con=np.log2(d[[c for c in cols if c.startswith(ln+'_Con')]]+1).mean(axis=1)
    for tag in ['6236','Taze','Combo']:
        t=np.log2(d[[c for c in cols if c.startswith(f'{ln}_{tag}')]]+1).mean(axis=1)
        lf=t-con; g=[s for s in SREBP if s in lf.index]
        print('  ',tag,'n=%d program %+.2f'%(len(g),lf.reindex(g).mean()),'HMGCS1 %+.2f HMGCR %+.2f MVD %+.2f'%(lf['HMGCS1'],lf['HMGCR'],lf['MVD']))
