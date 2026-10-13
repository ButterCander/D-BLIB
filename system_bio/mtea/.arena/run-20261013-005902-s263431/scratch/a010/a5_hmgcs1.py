import pandas as pd, numpy as np
B="D:/Cander/School/3-2/system_bio/mtea/"
conds = {"PI3Ki":"PI103","MEKi":"PD0325901","TAKi":"5Z","PI3Ki+TAKi":"PI103+5Z","PI3Ki+MEKi":"PI103+PD"}
genes = ["HMGCS1","HMGCR","MVD","MVK","SCD","SREBF1","FDPS","IDI1","INSIG1","LDLR","SQLE","ASNS","PHGDH","PSAT1","SHMT2","ODC1"]
rows=[]
for c,suf in conds.items():
    d = pd.read_csv(B+f"out/de_{suf}.tsv", sep="\t", index_col=0)
    for ge in genes:
        if ge in d.index:
            r=d.loc[ge]; rows.append((c,ge,r.log2FoldChange,r.lfc_shrunk,r.padj))
t = pd.DataFrame(rows, columns=["cond","gene","lfc_unshr","lfc_shr","padj"])
print("UNSHRUNK / SHRUNK (* padj<.05)")
for ge in genes:
    s=t[t.gene==ge].set_index("cond")
    print(f"{ge:7s}", " | ".join(f"{c}: {s.loc[c,'lfc_unshr']:+.2f}/{s.loc[c,'lfc_shr']:+.2f}{'*' if s.loc[c,'padj']<0.05 else ''}" for c in conds if c in s.index))
# authors' DEG lists (apeglm, |LFC|>=1 & padj<.05) from MOESM2
x = pd.ExcelFile(B+"si/MOESM2.xlsx")
print("\nauthors MOESM2 DEG lists (|LFC|>=1, padj<.05 only):")
for sh in x.sheet_names:
    d = x.parse(sh)
    d = d[[c for c in d.columns if not str(c).startswith("Unnamed")]]
    print(sh, len(d), "up", (d.log2FoldChange>0).sum(), "down", (d.log2FoldChange<0).sum(), "minabs", d.log2FoldChange.abs().min().round(3), "max padj", d.padj.max())
    for ge in ["HMGCS1","HMGCR","MVD","MVK","SCD","ASNS"]:
        r=d[d.geneSymbol==ge]
        if len(r): print("    ",ge, r.log2FoldChange.round(2).tolist(), r.padj.round(4).tolist())
