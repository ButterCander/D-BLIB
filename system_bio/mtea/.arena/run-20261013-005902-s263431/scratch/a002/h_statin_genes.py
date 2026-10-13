import pandas as pd
sym=pd.read_csv("pdac/statin_symbols.tsv",sep="\t").set_index("gene_id")["gene_short_name"]
for nm in ["PANC1_Simvastatin_High","MiaPaca2_Simvastatin_High","PANC1_Simvastatin_Low","MiaPaca2_Simvastatin_Low"]:
    d=pd.read_csv(f"out/statin_de_{nm}.tsv",sep="\t"); d["sym"]=d.gene.map(sym); d=d.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
    print(nm, " ".join(f"{g}:{d.loc[g,'logFC']:+.2f}{'*' if d.loc[g,'adj.P.Val']<.05 else ''}" for g in ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","SQLE","LDLR","INSIG1","SCD"] if g in d.index))
