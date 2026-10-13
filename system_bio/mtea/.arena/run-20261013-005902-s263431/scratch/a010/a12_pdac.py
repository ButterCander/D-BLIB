import pandas as pd, numpy as np, gzip
B="D:/Cander/School/3-2/system_bio/mtea/pdac/"
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TRUNK = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
out={}
for ln in ["CAPAN2","CFPAC1","MiaPaca2","PANC1"]:
    d = pd.read_csv(B+f"GSE238198_Pdac.{ln}.genes.normalizedCounts.txt.gz", sep="\t", index_col=0)
    cols = d.columns.tolist()
    print(ln, d.shape, cols)
    y = np.log2(d+1)
    grp = {}
    for c in cols:
        t = c.split("_")[1]; grp.setdefault(t, []).append(c)
    keep = (y.mean(axis=1) > np.log2(1+5))
    res={}
    for t in ["3h","12h","24h"]:
        lfc = y[grp[t]].mean(axis=1) - y[grp["NT"]].mean(axis=1)
        lfc = lfc[keep]
        res[t] = (lfc.reindex([g for g in SREBP if g in lfc.index]).mean(), lfc.reindex([g for g in TRUNK if g in lfc.index]).mean(),
                  lfc.reindex(["HMGCS1"]).iloc[0], lfc.reindex(["HMGCR"]).iloc[0], lfc.reindex(["MVK"]).iloc[0], lfc.reindex(["MVD"]).iloc[0], lfc.reindex(["SCD"]).iloc[0])
    print("  n by group", {k:len(v) for k,v in grp.items()}, "genes kept", int(keep.sum()))
    for t,v in res.items():
        print("   ",t,"program %.2f trunk %.2f | HMGCS1 %.2f HMGCR %.2f MVK %.2f MVD %.2f SCD %.2f"%v)
    # robustness: ratio of group-mean normalised counts with pseudo 1 (no log-then-mean)
    t="24h"
    r = np.log2((d[grp[t]].mean(axis=1)+1)/(d[grp["NT"]].mean(axis=1)+1))
    print("    alt (log ratio of mean counts) 24h: program %.2f"%r.reindex(SREBP).mean())
    out[ln]=res
