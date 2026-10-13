import pandas as pd, numpy as np, json
from scipy import stats
pd.set_option('display.width',250)
# AGS gene-level
files={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki+TAKi':'PI103+5Z','PI3Ki+MEKi':'PI103+PD'}
genes=['HMGCS1','HMGCR','MVD','MVK','SCD','SREBF1','FDPS','IDI1','INSIG1','LDLR','SQLE','ASNS']
rows=[]
for k,v in files.items():
    d=pd.read_csv(f'out/de_{v}.tsv',sep='\t').set_index('Symbol')
    for g in genes:
        if g in d.index:
            r=d.loc[g]; rows.append((k,g,round(r.log2FoldChange,2),round(r.lfc_shrunk,2),r.padj))
df=pd.DataFrame(rows,columns=['cond','gene','lfc','shrunk','padj'])
print(df.pivot(index='gene',columns='cond',values='lfc'))
print(df.pivot(index='gene',columns='cond',values='shrunk'))
print(df.pivot(index='gene',columns='cond',values='padj').map(lambda x:'*' if x<0.05 else 'ns'))
# DepMap pancreas
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0); cr=cr.loc[:,~cr.columns.duplicated()]
print('n pancreatic lines',cr.shape[0])
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','FNTA','FNTB','GGPS1','SQLE','LSS','DHCR7','FDFT1','IDI1','FDPS','SCAP','MBTPS1','RRM1','RRM2','CTPS1']:
    if g in cr: print(g, 'median %.2f'%cr[g].median(), 'frac<=-0.5 %.2f'%(cr[g]<=-0.5).mean(), 'n',cr[g].notna().sum())
cg=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0); cg=cg.loc[:,~cg.columns.duplicated()]
print('n gastric',cg.shape[0])
for g in ['HMGCS1','HMGCR','MVK','GGPS1','FNTA','SQLE','LSS']:
    print('gastric',g,'median %.2f'%cg[g].median())
