import pandas as pd, numpy as np
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0); cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas CRISPR rows (lines):',cr.shape[0])
md=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
print(md.loc[cr.index,['StrippedCellLineName','OncotreeSubtype']].head(50).to_string()[:3000])
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','FDFT1','SCD','SCAP','MBTPS1','SREBF2','RRM1','RRM2','CTPS1']:
    if g in cr.columns:
        s=cr[g].dropna(); print(f'{g:7} n={len(s)} median {s.median():+.2f}  frac<=-0.5 {(s<=-0.5).mean():.2f}')
    else: print(g,'missing')
cg=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0); cg=cg.loc[:,~cg.columns.duplicated()]
print('gastric lines',cg.shape[0])
for g in ['HMGCS1','HMGCR','MVK','SQLE','LSS']:
    print(g,'gastric median',round(cg[g].median(),2))
# AGS
ags=[i for i in cg.index if i=='ACH-000880']; print(ags)
if ags: print({g:round(cg.loc['ACH-000880',g],2) for g in ['HMGCS1','HMGCR','MVK','SQLE']})
