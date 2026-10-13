import pandas as pd, numpy as np
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
print(cr.shape); print(list(cr.index)[:5]); 
cr=cr.loc[:,~cr.columns.duplicated()]
print('n lines',len(cr))
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','FDFT1','SCD','SCAP','MBTPS1','SREBF2','MBTPS2','RRM1','RRM2','CTPS1','ASNS','PHGDH','ODC1','SRM','HMGCS2','GGPS1','RABGGTA','RABGGTB','PGGT1B','KRAS']:
    if g in cr.columns:
        v=cr[g].dropna(); print(f'{g:8} n={len(v)} median={v.median():.2f} frac<=-0.5={(v<=-0.5).mean():.2f} min={v.min():.2f} max={v.max():.2f}')
    else: print(g,'missing')
gs=pd.read_csv('external/crispr_gene_summary.tsv',sep='\t'); print(gs.head()); print(gs.shape)
gr=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0); gr=gr.loc[:,~gr.columns.duplicated()]
print('gastric lines',len(gr))
for g in ['HMGCS1','HMGCR','MVK','GGPS1','FNTA','SQLE','LSS']:
    print(g,'gastric median',round(gr[g].median(),2))
