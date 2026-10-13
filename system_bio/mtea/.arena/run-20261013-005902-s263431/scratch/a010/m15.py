import sys, pandas as pd, numpy as np
from scipy.stats import binomtest
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
for line in ['PSN1','SW1990']:
    d=pd.read_csv(base+f'pdac/GSE315130_{line}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t')
    d=d.drop_duplicates('Symbol').set_index('Symbol')
    for cond in ['6236','Taze','Combo']:
        lf=np.log2(d[[f'{line}_{cond}_{i}' for i in (1,2,3)]]+1).mean(axis=1)-np.log2(d[[f'{line}_Con_{i}' for i in (1,2,3)]]+1).mean(axis=1)
        pr=lf.reindex([g for g in SREBP if g in lf.index]).dropna()
        sp=binomtest(int((pr<0).sum()),len(pr),0.5).pvalue
        print(line,cond,'program mean %+.2f n=%d down=%d sign p=%.3f'%(pr.mean(),len(pr),(pr<0).sum(),sp),{g:round(lf.get(g,np.nan),2) for g in ['HMGCS1','HMGCR','MVK','MVD']})
