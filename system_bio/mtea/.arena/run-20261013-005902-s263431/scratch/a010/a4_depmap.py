import pandas as pd, numpy as np
B="D:/Cander/School/3-2/system_bio/mtea/external/"
s = pd.read_csv(B+"crispr_gene_summary.tsv", sep="\t", index_col=0)
s.index = [i.split(" (")[0] for i in s.index]
c = pd.read_csv(B+"crispr_pancreas_rows.tsv", sep="\t", index_col=0)
c.columns = [x.split(" (")[0] for x in c.columns]
c = c.loc[:, ~c.columns.duplicated()]
g = pd.read_csv(B+"crispr_gastric_rows.tsv", sep="\t", index_col=0)
g.columns = [x.split(" (")[0] for x in g.columns]
g = g.loc[:, ~g.columns.duplicated()]
genes = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","GGPS1","FNTA","FNTB","SQLE","LSS","DHCR7","FDFT1","CYP51A1","SCD","SCAP","SREBF2","MBTPS1","RRM1","RRM2","CTPS1","ASNS","PSAT1","PHGDH","ODC1","COQ3","TXNRD1"]
print("gene  panc_med panc_frac(<=-0.5)  gastr_med | pan(all 1178 lines) mean  frac_dep")
for ge in genes:
    if ge in c.columns:
        print(f"{ge:8s} {c[ge].median():6.2f} {(c[ge]<=-0.5).mean():5.2f} (n={c[ge].notna().sum()}) | {g[ge].median():6.2f} | {s.loc[ge,'mean']:6.2f} {s.loc[ge,'frac_dep']:5.2f}")
print("n panc rows", c.shape[0], "gastric", g.shape[0])
# what fraction of all genes have frac_dep >= .9 (common essential)? 
print("genes with pan frac_dep>=0.9:", (s.frac_dep>=0.9).sum(), "of", len(s))
# AGS
ags = "ACH-000880"
print("AGS in gastric rows?", ags in g.index)
for ge in ["HMGCS1","HMGCR","MVK","MVD","SQLE","ASNS","PHGDH","ODC1","RRM1","RRM2","COQ3"]:
    print(ge, g.loc[ags,ge] if ags in g.index else None)
