import pandas as pd
d=pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv")
gi=pd.read_csv("tsirvouli/gene_info.tsv",sep="\t",index_col=0)
d['sym']=d.Gene_name.map(gi.gene_name)
genes=["HMGCS1","HMGCR","MVD","MVK","IDI1","FDPS","SCD","SREBF1","SREBF2","INSIG1","SCAP","LDLR","SQLE"]
for tp in ['8h','24h']:
    print(tp)
    for g in genes:
        row=[]
        for c in ['PI','PD','OXO','PIPD','PIOXO']:
            r=d[(d.sym==g)&(d.Contrast==f'Diff_{c}_{tp}')]
            row.append(f"{r.logFC.iloc[0]:+.2f}{'*' if r['adj.P.Val'].iloc[0]<.05 else ' '}" if len(r) else 'NA')
        print(f"  {g:7}", *row)
