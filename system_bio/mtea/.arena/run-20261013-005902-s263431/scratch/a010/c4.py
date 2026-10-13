import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
names={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki+TAKi':'PI103+5Z','PI3Ki+MEKi':'PI103+PD'}
de={k:pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0) for k,f in names.items()}
genes=['HMGCS1','HMGCR','MVD','MVK','SCD','SREBF1','FDPS','IDI1','ASNS','SQLE','LDLR','INSIG1']
rows=[]
for g in genes:
    r=[g]
    for k,d in de.items():
        if g in d.index:
            x=d.loc[g]; r.append(f"{x.lfc_shrunk:+.2f}{'*' if x.padj<0.05 else ''}/{x.log2FoldChange:+.2f}")
        else: r.append('NA')
    rows.append(r)
print(pd.DataFrame(rows,columns=['gene']+list(names)).to_string(index=False))
