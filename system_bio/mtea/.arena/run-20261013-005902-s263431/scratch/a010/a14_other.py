import numpy as np, pandas as pd, gzip
B="D:/Cander/School/3-2/system_bio/mtea/"
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
# ---- GSE206267 raw counts
c = pd.read_csv(B+"pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz", sep="\t", index_col=0)
print("GSE206267", c.shape, c.columns.tolist())
gi = pd.read_csv(B+"tsirvouli/gene_info.tsv", sep="\t", index_col=0)["gene_name"]
c = c[c.sum(axis=1)>=30]
c["sym"] = gi.reindex(c.index).values
c = c.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
nc = c / c.sum() * c.sum().median(); L = np.log2(nc+1)
for line,pre in (("KP4","KP4"),("MiaPaCa2","Mia")):
    d = L[[f"{pre}-T100-{i}" for i in (1,2,3)]].mean(axis=1) - L[[f"{pre}-DMSO-{i}" for i in (1,2,3)]].mean(axis=1)
    pr = d.reindex([g for g in SREBP if g in d.index])
    print(f"  {line}: program mean {pr.mean():+.2f} (n={pr.notna().sum()}), up {(pr>0).sum()}/{pr.notna().sum()}; HMGCS1 {d['HMGCS1']:+.2f} HMGCR {d['HMGCR']:+.2f} MVD {d['MVD']:+.2f} MVK {d['MVK']:+.2f}")
# ---- GSE279749 FPKM
d = pd.read_csv(B+"pdac/GSE279749_gene_fpkm_cell_lines.txt.gz", sep="\t")
print("GSE279749", d.shape, [x for x in d.columns][:40])
d = d.drop_duplicates("gene_name").set_index("gene_name")
cols = d.columns.tolist()
print(cols)
