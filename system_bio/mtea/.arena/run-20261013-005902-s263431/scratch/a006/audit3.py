import pandas as pd, numpy as np, json, gzip
B='D:/Cander/School/3-2/system_bio/mtea/'
cr=pd.read_csv(B+'external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas rows',cr.shape)
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','SCD','RRM1','RRM2','CTPS1','SCAP','MBTPS1','SREBF2','INSIG1','LDLR','HMGCS2']:
    if g in cr.columns:
        v=cr[g].dropna()
        print('%-7s n=%d median=%.2f  frac<=-0.5=%.2f  frac<=-1=%.2f min=%.2f max=%.2f'%(g,len(v),v.median(),(v<=-0.5).mean(),(v<=-1).mean(),v.min(),v.max()))
    else: print(g,'missing')
# gastric rows
cg=pd.read_csv(B+'external/crispr_gastric_rows.tsv',sep='\t',index_col=0,nrows=None)
cg=cg.loc[:,~cg.columns.duplicated()]
print('gastric rows',cg.shape)
for g in ['HMGCS1','HMGCR','MVK','SQLE','GGPS1','FNTA']:
    print(g,'AGS',cg.loc['ACH-000880',g] if 'ACH-000880' in cg.index else 'noAGS','gastric median %.2f'%cg[g].median())
# how many pancreas lines (nonNA) for HMGCS1 and pan-essentiality across all DepMap (gene summary)
gs=pd.read_csv(B+'external/crispr_gene_summary.tsv',sep='\t',nrows=3)
print(gs.head(3).T.head(12))
