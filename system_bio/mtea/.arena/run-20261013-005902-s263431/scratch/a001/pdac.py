import numpy as np, pandas as pd, sys
sys.stdout.reconfigure(encoding="utf-8")
L = ["CAPAN2","CFPAC1","MiaPaca2","PANC1"]
G = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","SCD","SQLE","LSS","INSIG1","LDLR","SREBF2","RRM1","RRM2","CTPS1","TYMS","FNTA","GGPS1"]
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
d0 = pd.read_csv("out/pdac_de_CAPAN2_h24.tsv", sep="\t"); print(d0.columns.tolist(), len(d0))
lf = {}; pa = {}; allgenes={}
for l in L:
    d = pd.read_csv(f"out/pdac_de_{l}_h24.tsv", sep="\t").set_index("gene")
    lf[l]=d["logFC"]; pa[l]=d["adj.P.Val"]; allgenes[l]=d
LF=pd.DataFrame(lf); PA=pd.DataFrame(pa)
print("== log2FC trametinib 24h (padj<.05 marked *)")
t = LF.reindex(G).round(2).astype(str) + np.where(PA.reindex(G)<0.05,"*"," ")
print(t.to_string())
print("program mean (24 SREBP genes):", {l: round(LF.reindex(SREBP)[l].mean(),2) for l in L})
print("n genes up of those measured:", {l:(int((LF.reindex(SREBP)[l]>0).sum()),int(LF.reindex(SREBP)[l].notna().sum())) for l in L})
# rule: >=0.5 and padj<.05 in >=3 of 4
up = ((LF>=0.5)&(PA<0.05)).sum(axis=1)
print("genes with >=0.5 & padj<.05 in >=3 lines:", int((up>=3).sum()))
print(up.reindex(["HMGCS1","HMGCR","MVK","MVD","SCD","IDI1","FDPS"]).to_dict())
print("n DEGs padj<.05:", {l:int((PA[l]<0.05).sum()) for l in L})
# per-line number of key genes >= .5
for g in ["HMGCS1","HMGCR","MVK","MVD","SCD"]:
    print(g, [ (l, round(LF.loc[g,l],2), round(PA.loc[g,l],4)) for l in L])
