import numpy as np, pandas as pd
B="D:/Cander/School/3-2/system_bio/mtea/"
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
d = pd.read_csv(B+"pdac/GSE279749_gene_fpkm_cell_lines.txt.gz", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
for nm,k in (("AsPC-1","A"),("PANC-1","P"),("CFPAC-1","C"),("SW1990","S9")):
    lf = np.log2(d[[f"{k}_M_{i}" for i in (1,2,3)]]+1).mean(axis=1) - np.log2(d[[f"{k}_NT_{i}" for i in (1,2,3)]]+1).mean(axis=1)
    pr = lf.reindex([g for g in SREBP if g in lf.index]).dropna()
    print(nm, "program %+.2f (n=%d, up %d)"%(pr.mean(), len(pr), (pr>0).sum()), "HMGCS1 %+.2f HMGCR %+.2f MVD %+.2f"%(lf.get("HMGCS1",np.nan), lf.get("HMGCR",np.nan), lf.get("MVD",np.nan)))
# how many genes; MAPK-pathway positive control? DUSP6/ETV4/SPRY pathway genes (do they fall with MRTX1133?)
for nm,k in (("AsPC-1","A"),("PANC-1","P"),("CFPAC-1","C"),("SW1990","S9")):
    lf = np.log2(d[[f"{k}_M_{i}" for i in (1,2,3)]]+1).mean(axis=1) - np.log2(d[[f"{k}_NT_{i}" for i in (1,2,3)]]+1).mean(axis=1)
    print(nm, {g: round(float(lf.get(g,np.nan)),2) for g in ["DUSP6","ETV4","ETV5","SPRY4","FOSL1","MYC","CCND1","DUSP4"]})
