import sys, pandas as pd, numpy as np
from scipy import stats
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
c=pd.read_csv(base+'pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
print(c.shape, list(c.columns))
c=c[c.sum(axis=1)>=30]
gi=pd.read_csv(base+'tsirvouli/gene_info.tsv',sep='\t',index_col=0)['gene_name']
c['sym']=gi.reindex(c.index).values
c=c.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym')
nc=c/c.sum()*c.sum().median(); L=np.log2(nc+1)
for line,pre in (('KP4','KP4'),('MiaPaCa2','Mia')):
    cols_t=[f'{pre}-T100-{i}' for i in (1,2,3)]; cols_c=[f'{pre}-DMSO-{i}' for i in (1,2,3)]
    d=L[cols_t].mean(axis=1)-L[cols_c].mean(axis=1)
    # welch t per gene for HMGCS1 etc
    out={}
    for g in ['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SQLE']:
        t,p=stats.ttest_ind(L.loc[g,cols_t],L.loc[g,cols_c],equal_var=False)
        out[g]=f'{d[g]:+.2f}(p={p:.3f})'
    prog=d.reindex([g for g in SREBP if g in d.index])
    print(line,'program mean %.2f n=%d up %d'%(prog.mean(),prog.notna().sum(),(prog>0).sum()),out)
    from scipy.stats import binomtest
    print('  sign test p',binomtest(int((prog>0).sum()),int(prog.notna().sum()),0.5).pvalue)
