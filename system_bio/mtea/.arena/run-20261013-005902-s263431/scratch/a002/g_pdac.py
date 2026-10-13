import numpy as np, pandas as pd
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TRUNK=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
LINES=["CAPAN2","CFPAC1","MiaPaca2","PANC1"]
D={(l,t):pd.read_csv(f"out/pdac_de_{l}_{t}.tsv",sep="\t").set_index("gene") for l in LINES for t in ("h3","h12","h24")}
print("24h logFC (padj<.05 *)")
for g in ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","SCD","SQLE","LDLR","INSIG1","SREBF2","SREBF1"]:
    row=[]
    for l in LINES:
        d=D[(l,"h24")]
        if g in d.index:
            r=d.loc[g]; row.append(f"{r.logFC:+.2f}{'*' if r['adj.P.Val']<.05 else ' '}")
        else: row.append("NA")
    print(f"  {g:7}",*row)
print("program mean (24 genes) and #up, trunk(6) mean, by time")
for l in LINES:
    out=[]
    for t in ("h3","h12","h24"):
        d=D[(l,t)]["logFC"]
        x=d.reindex(SREBP).dropna()
        out.append(f"{t}: prog {x.mean():+.2f} ({int((x>0).sum())}/{len(x)} up) trunk {d.reindex(TRUNK).mean():+.2f}")
    print(f"  {l:9}"," | ".join(out))
print("DEG counts padj<.05 at 24h:",{l:int((D[(l,'h24')]['adj.P.Val']<.05).sum()) for l in LINES})
print("n >=0.5 & padj<.05 in lines (24h) for trunk genes:")
for g in TRUNK+["SCD"]:
    n=sum(1 for l in LINES if g in D[(l,'h24')].index and D[(l,'h24')].loc[g,'logFC']>=0.5 and D[(l,'h24')].loc[g,'adj.P.Val']<.05)
    print(" ",g,n)
