import pandas as pd
conds=[("PI3Ki","PI103"),("MEKi","PD0325901"),("TAKi","5Z"),("PI3Ki+TAKi","PI103+5Z"),("PI3Ki+MEKi","PI103+PD")]
genes=["HMGCS1","HMGCR","MVD","MVK","IDI1","FDPS","SCD","SREBF1","SREBF2","SQLE","LDLR","INSIG1","TDO2","ASNS"]
rows={}
for n,c in conds:
    d=pd.read_csv(f"out/de_{c}.tsv",sep="\t",index_col=0)
    d=d[~d.index.duplicated()]
    rows[n]=d
print("shrunk (lfc_shrunk) | unshrunk (log2FoldChange); * padj<.05")
for g in genes:
    s=[]
    for n,_ in conds:
        d=rows[n]
        if g in d.index:
            r=d.loc[g]; s.append(f"{r.lfc_shrunk:+.2f}/{r.log2FoldChange:+.2f}{'*' if r.padj<.05 else ' '}")
        else: s.append("NA")
    print(f"{g:7}", " | ".join(s))
