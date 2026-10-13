import pandas as pd, numpy as np
pd.set_option('display.width',250)
d=pd.read_csv('out/pdac_r2d.tsv',sep='\t',index_col=0)
iso=["HMGCR","HMGCS1","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","DHCR7","DHCR24","SC5D","LDLR","INSIG1","ACAT2","NSDHL","MSMO1","HSD17B7","EBP","TM7SF2","SREBF2","SCD"]
print(d.reindex([g for g in iso if g in d.index])[['mean_lfc','n_lines','n_up','pdac_median','pdac_frac_dep','pan_mean','ags_gse','ags_tsir','statin','class']].round(2).to_string())
print(d['class'].value_counts())
# per-line values for HMGCR etc at 24h
for l in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    x=pd.read_csv(f'out/pdac_de_{l}_h24.tsv',sep='\t').set_index('gene')
    print(l, {g:(round(x.loc[g,'logFC'],2),round(x.loc[g,'adj.P.Val'],3)) for g in ['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','LDLR','INSIG1','SCD','SQLE'] if g in x.index})
