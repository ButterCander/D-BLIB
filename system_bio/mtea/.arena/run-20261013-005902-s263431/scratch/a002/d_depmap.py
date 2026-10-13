import pandas as pd, numpy as np
cr=pd.read_csv("external/crispr_pancreas_rows.tsv",sep="\t",index_col=0)
cr=cr.loc[:,~cr.columns.duplicated()]
print(cr.shape)
print(cr.index[:5].tolist())
mod=pd.read_csv("external/Model_24Q4.csv")
print(mod.columns.tolist()[:30])
for g in ["HMGCS1","HMGCR","MVK","MVD","PMVK","IDI1","FDPS","GGPS1","FNTA","FNTB","SQLE","LSS","DHCR7","SCD","SREBF2","SCAP","MBTPS1","MBTPS2","INSIG1","LDLR","RRM1","RRM2","CTPS1"]:
    if g in cr.columns:
        v=cr[g].dropna()
        print(f"{g:7} n={len(v)} median {v.median():+.2f}  frac<=-0.5 {(v<=-0.5).mean():.2f}  frac<=-1 {(v<=-1).mean():.2f}")
    else: print(g,'absent')
