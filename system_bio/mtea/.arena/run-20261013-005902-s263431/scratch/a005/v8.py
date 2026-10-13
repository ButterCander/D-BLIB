import pandas as pd, numpy as np, gzip
from scipy import stats
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
d=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t')
print(d.shape, d.columns.tolist()[:14])
for f in ['GSE315130_PSN1_cell_line_gene_sample_FPKM_with_symbol.txt.gz','GSE315130_SW1990_cell_line_gene_sample_FPKM_with_symbol.txt.gz']:
    x=pd.read_csv('pdac/'+f,sep='\t'); print(f, x.shape, x.columns.tolist()[:20])
