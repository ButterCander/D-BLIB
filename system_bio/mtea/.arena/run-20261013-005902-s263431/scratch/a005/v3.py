import pandas as pd, numpy as np, json, re
from scipy import stats
pd.set_option('display.width',250)
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv'); gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
print(sorted(off.Contrast.unique())[:40])
tr=['HMGCS1','HMGCR','MVK','MVD','IDI1','INSIG1','LDLR','SCD','SQLE']
for t in ['1h','2h','4h','8h','24h']:
    r={}
    for c in ['PI','PD','PIPD','OXO','PIOXO']:
        x=off[(off.Contrast==f'Diff_{c}_{t}')&off.sym.isin(tr)].set_index('sym')
        r[c]=x.logFC.reindex(tr).round(2)
    print(t); print(pd.DataFrame(r))
# full-dose design: were these p adj?
x=off[(off.Contrast=='Diff_PIPD_24h')&off.sym.isin(['HMGCS1'])]; print(x[['logFC','adj.P.Val']])
x=off[(off.Contrast=='Diff_PI_24h')&off.sym.isin(['HMGCS1','MVD','IDI1'])]; print(x[['sym','logFC','adj.P.Val']])
