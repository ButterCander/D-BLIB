import pandas as pd, numpy as np
conds={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki+TAKi':'PI103+5Z','PI3Ki+MEKi':'PI103+PD'}
genes=['HMGCS1','HMGCR','MVD','MVK','FDPS','IDI1','SCD','SREBF1','SREBF2','INSIG1','LDLR','SQLE','LSS','DHCR7','ASNS']
rows=[]
for k,f in conds.items():
    d=pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0)
    for g in genes:
        if g in d.index:
            r=d.loc[g]; rows.append((g,k,round(r.log2FoldChange,2),round(r.lfc_shrunk,2),'*' if r.padj<0.05 else '', round(r.padj,4) if pd.notna(r.padj) else None, round(r.baseMean,0)))
df=pd.DataFrame(rows,columns=['gene','cond','mle','shrunk','sig','padj','baseMean'])
print(df.pivot(index='gene',columns='cond',values='shrunk').round(2)[list(conds)])
print(df.pivot(index='gene',columns='cond',values='mle').round(2)[list(conds)])
print(df.pivot(index='gene',columns='cond',values='sig')[list(conds)])
print(df.pivot(index='gene',columns='cond',values='padj')[list(conds)])
