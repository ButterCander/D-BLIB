import pandas as pd, numpy as np
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas lines',cr.shape)
print(cr.iloc[:3,:5])
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','FDFT1','SCAP','MBTPS1','SREBF2','SREBF1','SCD','FASN','LDLR','INSIG1','RRM1','RRM2','CTPS1','GMPS','PPAT']:
    if g in cr.columns:
        s=cr[g].dropna()
        print('%-7s n=%d median %.2f  frac<=-0.5 %.2f  (n dep %d)  min %.2f max %.2f'%(g,len(s),s.median(),(s<=-0.5).mean(),(s<=-0.5).sum(),s.min(),s.max()))
    else: print(g,'missing')
gs=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0); gs=gs.loc[:,~gs.columns.duplicated()]
print('gastric lines',gs.shape)
for g in ['HMGCS1','HMGCR','MVK','GGPS1','FNTA','SQLE','LSS']:
    print(g,'gastric median %.2f'%gs[g].median())
print(gs.index[:5].tolist())
