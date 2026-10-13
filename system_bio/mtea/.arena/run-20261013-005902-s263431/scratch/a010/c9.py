import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
print(d.shape)
for k,nm in (('A','AsPC-1'),('P','PANC-1'),('C','CFPAC-1'),('S9','SW1990')):
    a=np.log2(d[[f'{k}_M_{i}' for i in (1,2,3)]]+1).mean(axis=1); b=np.log2(d[[f'{k}_NT_{i}' for i in (1,2,3)]]+1).mean(axis=1)
    lf=a-b; s=lf.reindex([g for g in SREBP if g in lf.index]).dropna()
    rng=np.random.default_rng(1); pool=lf.dropna().values
    rnd=np.array([rng.choice(pool,len(s),replace=False).mean() for _ in range(3000)])
    t=lf.reindex(TR)
    print(nm,len(s),int((s>0).sum()),round(s.mean(),2),round((rnd>=s.mean()).mean(),3),' '.join(f"{g}:{t[g]:+.2f}" for g in TR if pd.notna(t[g])), ' genes:',len(lf.dropna()))
    # background baseline: what fraction of low-count genes are retained
print(d['gene_name'].notna().sum() if 'gene_name' in d else '')
print(d.filter(like='_NT_').mean(axis=1).describe())
