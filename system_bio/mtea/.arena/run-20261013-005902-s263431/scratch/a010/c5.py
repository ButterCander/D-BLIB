import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas rows',cr.shape)
cg=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0)
cg.columns=[c.split(' (')[0] for c in cg.columns]; cg=cg.loc[:,~cg.columns.duplicated()]
print('gastric rows',cg.shape)
summ=pd.read_csv('external/crispr_gene_summary.tsv',sep='\t',index_col=0); summ.index=[c.split(' (')[0] for c in summ.index]
genes=['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','FNTA','FNTB','GGPS1','SQLE','LSS','DHCR7','FDFT1','SCD','RRM1','RRM2','CTPS1','SCAP','SREBF2','LDLR','ASNS']
rows=[]
for g in genes:
    if g not in cr.columns: rows.append((g,'NA')); continue
    rows.append((g, round(cr[g].median(),2), f"{(cr[g]<=-0.5).sum()}/{cr[g].notna().sum()}", round((cr[g]<=-0.5).mean(),3), round(cg[g].median(),2) if g in cg.columns else None, round(summ.loc[g,'mean'],2), round(summ.loc[g,'frac_dep'],2)))
print(pd.DataFrame(rows,columns=['gene','panc_med','n_dep','frac','gastric_med','pan_mean(1178)','pan_fracdep']).to_string(index=False))
