import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")
conds = {"PI3Ki":"PI103","MEKi":"PD0325901","TAKi":"5Z","PI3Ki+TAKi":"PI103+5Z","PI3Ki+MEKi":"PI103+PD"}
G = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","SCD","SREBF1","SREBF2","SQLE","LSS","INSIG1","LDLR","ASNS","TDO2","ODC1","COQ3","TXNRD1"]
for col in ["log2FoldChange","lfc_shrunk"]:
    print("==",col)
    rows=[]
    for k,v in conds.items():
        d=pd.read_csv(f"out/de_{v}.tsv",sep="\t",index_col=0)
        rows.append(d.reindex(G)[col].rename(k))
    print(pd.concat(rows,axis=1).round(2).to_string())
print("== padj")
rows=[]
for k,v in conds.items():
    d=pd.read_csv(f"out/de_{v}.tsv",sep="\t",index_col=0)
    rows.append(d.reindex(G)["padj"].rename(k))
print(pd.concat(rows,axis=1).round(4).to_string())
d=pd.read_csv("out/de_PI103.tsv",sep="\t",index_col=0); print(d.columns.tolist(), len(d))
