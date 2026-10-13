import pandas as pd, numpy as np, json
from scipy import stats
print("=== A. GSE285616 raw-count log2FC (batch-paired, my own normalisation) ===")
c=pd.read_csv('GSE285616_counts.csv.gz',index_col=0); c=c[~c.index.duplicated()]
x=c.values.astype(float); pos=(x>0).all(axis=1)
sf=np.median(x[pos]/np.exp(np.log(x[pos]).mean(axis=1,keepdims=True)),axis=0)
n=pd.DataFrame(x/sf,index=c.index,columns=c.columns)
L=np.log2(n+1)
order=['ctrl','PI','PD','5Z','PI5Z','PIPD']
def lfc(g,k):
    # mean over batches of log2((norm+1)) difference vs control of same batch
    v=[L.loc[g].iloc[b*6+k]-L.loc[g].iloc[b*6] for b in range(4)]
    return np.mean(v)
def lfc2(g,k):
    a=n.loc[g].iloc[[b*6+k for b in range(4)]].mean(); b0=n.loc[g].iloc[[b*6 for b in range(4)]].mean()
    return np.log2(a/b0)
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','SREBF1','ASNS']:
    print(g,{order[k]:round(lfc2(g,k),2) for k in range(1,6)})
print("de_*.tsv (apeglm shrunk / unshrunk / padj)")
for f,lab in [('PI103','PI'),('PD0325901','PD'),('5Z','5Z'),('PI103+5Z','PI5Z'),('PI103+PD','PIPD')]:
    d=pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0)
    print(lab,{g:(round(d.loc[g,'lfc_shrunk'],2),round(d.loc[g,'log2FoldChange'],2),'*' if d.loc[g,'padj']<0.05 else 'ns') for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','SREBF1']})
print()
print("=== C. DepMap ===")
cr=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
cr.columns=[c_.split(' (')[0] for c_ in cr.columns]; cr=cr.loc[:,~cr.columns.duplicated()]
print('pancreas rows',cr.shape)
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','FNTA','FNTB','GGPS1','SQLE','LSS','DHCR7','SCD','RRM1','RRM2']:
    if g in cr: print(g,'median',round(cr[g].median(),2),'frac<=-0.5',round((cr[g]<=-0.5).mean(),2))
gr=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0)
gr.columns=[c_.split(' (')[0] for c_ in gr.columns]; gr=gr.loc[:,~gr.columns.duplicated()]
print('gastric rows',gr.shape,'AGS' , 'ACH-000880' in gr.index)
for g in ['HMGCS1','HMGCR','MVK','SQLE']:
    print(g,'gastric median',round(gr[g].median(),2),'AGS',round(gr.loc['ACH-000880',g],2) if 'ACH-000880' in gr.index else None)
print()
print("=== D. Tsirvouli official limma ===")
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv'); gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
for g in ['HMGCS1','HMGCR','MVK','MVD']:
    out={}
    for c_ in ['PI','PD','OXO','PIPD','PIOXO']:
        for h in ['4h','8h','24h']:
            s=off[(off.Contrast==f'Diff_{c_}_{h}')&(off.sym==g)]
            if len(s): out[f'{c_}_{h}']=(round(s.logFC.iloc[0],2),'*' if s['adj.P.Val'].iloc[0]<0.05 else 'ns')
    print(g,out)
print(sorted(off.Contrast.unique())[:40])
