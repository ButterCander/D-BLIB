import pandas as pd, numpy as np, json
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0); cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas lines',cr.shape)
for g in ['HMGCS1','HMGCR','MVK','MVD','PMVK','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','SCD','SCAP','SREBF2','MBTPS1','LDLR','RRM1','RRM2','CTPS1']:
    if g in cr.columns: print(g,'median %.2f'%cr[g].median(),'frac<=-0.5 %.2f'%((cr[g]<=-0.5).mean()), 'n',cr[g].notna().sum())
    else: print(g,'missing')
gs=pd.read_csv('external/crispr_gene_summary.tsv',sep='\t'); print(gs.shape); print(gs.head(3).to_string()); 
print(gs[gs.iloc[:,0].isin(['HMGCS1','HMGCR','SQLE','MVK','MVD','GAPDH','RPL3'])].to_string())
gr=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0); gr=gr.loc[:,~gr.columns.duplicated()]; print('gastric',gr.shape)
for g in ['HMGCS1','HMGCR','SQLE']: print(g,'gastric median %.2f'%gr[g].median(),'frac %.2f'%((gr[g]<=-0.5).mean()))
