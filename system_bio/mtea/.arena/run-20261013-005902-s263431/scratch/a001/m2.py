import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding="utf-8")
x = pd.ExcelFile("si/MOESM2.xlsx")
A={}
for s in x.sheet_names:
    d=pd.read_excel(x,sheet_name=s); A[s]=d
    print(s,len(d),"min|lfc|=%.2f"%d.log2FoldChange.abs().min(), "max padj=%.3f"%d.padj.max(), "up",(d.log2FoldChange>0).sum(),"down",(d.log2FoldChange<0).sum())
G=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS","SQLE","LSS","SCD","INSIG1","LDLR","ASNS","TDO2","MVD"]
out=pd.DataFrame({s:A[s].drop_duplicates("geneSymbol").set_index("geneSymbol").reindex(G).log2FoldChange for s in A})
print(out.round(2).to_string())
# TDO2 id
for s in A:
    d=A[s]; print(s, d[d.geneSymbol=="TDO2"][["geneID","log2FoldChange","padj"]].values.tolist())
