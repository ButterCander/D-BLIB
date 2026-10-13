import pandas as pd, numpy as np
B="D:/Cander/School/3-2/system_bio/mtea/"
x = pd.ExcelFile(B+"si/MOESM2.xlsx")
d = x.parse("MEKi"); d = d[[c for c in d.columns if not str(c).startswith("Unnamed")]]
print(d[d.geneSymbol.isin(["TDO2","IDO1","IDO2"])])
hg = pd.read_csv(B+"mteapy/src/mteapy/data/HumanGEM_genes.tsv", sep="\t")
print(hg[hg.geneSymbol.isin(["TDO2","IDO1","IDO2"])])
# count symbols in the authors' DEG lists whose Ensembl ID differs from HumanGEM's for same symbol
hgm = hg.drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
tot=set()
for sh in x.sheet_names:
    dd = x.parse(sh); dd = dd[[c for c in dd.columns if not str(c).startswith("Unnamed")]]
    m = dd[dd.geneSymbol.isin(hgm.index)]
    bad = m[m.geneID != hgm.reindex(m.geneSymbol).values]
    tot |= set(bad.geneSymbol)
print("symbols with different ID than Human-GEM (union over 5 sheets):", len(tot))
# raw counts of TDO2 in GEO counts
c = pd.read_csv(B+"GSE285616_counts.csv.gz", index_col=0)
print(c.columns.tolist()[:30], c.shape)
print(c.loc[c.index.isin(["TDO2","IDO1"])].T if c.index.isin(["TDO2"]).any() else "TDO2 not index")
