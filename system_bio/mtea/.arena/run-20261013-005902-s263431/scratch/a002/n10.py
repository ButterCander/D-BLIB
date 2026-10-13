import numpy as np, pandas as pd
from scipy import stats
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
pm={}
for l in ("CAPAN2","CFPAC1","MiaPaca2","PANC1"):
    d=pd.read_csv(f'out/pdac_de_{l}_h24.tsv',sep='\t').set_index('gene')
    s=d.reindex([g for g in SREBP if g in d.index])
    pm[l]=s.logFC.mean(); 
    print(l,'program mean',round(s.logFC.mean(),2),'n',len(s),'up',(s.logFC>0).sum(),'trunk mean',round(d.reindex(TR).logFC.mean(),2),'4-gene (HMGCS1,HMGCR,MVK,MVD)',round(d.reindex(['HMGCS1','HMGCR','MVK','MVD']).logFC.mean(),2))
# MRTX1133
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
print([c for c in d.columns][:40])
