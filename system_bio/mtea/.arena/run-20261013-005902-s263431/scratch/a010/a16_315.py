import numpy as np, pandas as pd
from scipy import stats
B="D:/Cander/School/3-2/system_bio/mtea/pdac/"
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
for ln in ["PSN1","SW1990"]:
    d = pd.read_csv(B+f"GSE315130_{ln}_cell_line_gene_sample_FPKM_with_symbol.txt.gz", sep="\t")
    d = d.dropna(subset=["Symbol"]).drop_duplicates("Symbol").set_index("Symbol")
    d = d[[c for c in d.columns if c.startswith(ln)]]
    # filter low-expression
    keep = d.mean(axis=1) > 1
    for arm in ["6236","Combo","Taze"]:
        lf = np.log2(d[[f"{ln}_{arm}_{i}" for i in (1,2,3)]]+1).mean(axis=1) - np.log2(d[[f"{ln}_Con_{i}" for i in (1,2,3)]]+1).mean(axis=1)
        lf = lf[keep]
        pr = lf.reindex([g for g in SREBP if g in lf.index]).dropna()
        sg = stats.binomtest(int((pr<0).sum()), len(pr), 0.5).pvalue
        print(ln, arm, "program %+.2f n=%d down %d up %d (sign test p %.3f)"%(pr.mean(), len(pr), (pr<0).sum(), (pr>0).sum(), sg),
              "| DUSP6 %+.2f FOSL1 %+.2f MYC %+.2f HMGCS1 %+.2f"%(lf.get("DUSP6",np.nan), lf.get("FOSL1",np.nan), lf.get("MYC",np.nan), lf.get("HMGCS1",np.nan)))
