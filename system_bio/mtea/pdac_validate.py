"""Independent replication of the mevalonate-trunk induction and a program-level (multi-gene) test.
(1) GSE206267: KP4 and MiaPaCa-2, trametinib 100 nM(?) vs DMSO 48 h, raw counts (limma-trend on log2(CPM+1)... here log2 normalised counts + 1).
(2) GSE315130: PSN1 and SW1990, RAS(ON) inhibitor RMC-6236 (daraxonrasib class) vs DMSO, FPKM (log2(FPKM+1), Welch t per gene).
(3) Program score: SREBP2/sterol program (pre-listed 24 genes) across all datasets; sign counts and mean log2FC."""
import numpy as np, pandas as pd
from scipy import stats
TRUNK = ["HMGCS1", "HMGCR", "MVK", "MVD", "IDI1", "FDPS"]
SREBP = ["HMGCS1", "HMGCR", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1", "MSMO1", "NSDHL", "HSD17B7", "DHCR7", "DHCR24", "SC5D", "EBP", "TM7SF2", "ACAT2", "INSIG1", "LDLR", "SCD", "FASN"]
BRANCH = ["SQLE", "LSS", "CYP51A1", "MSMO1", "NSDHL", "HSD17B7", "DHCR7", "DHCR24", "SC5D", "EBP", "TM7SF2"]
res = {}
# ---- (1) GSE206267 raw counts: library-size normalisation, log2(x+1), per-gene Welch t (3 vs 3)
c = pd.read_csv("pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz", sep="\t", index_col=0)
hg = pd.read_csv("mteapy/src/mteapy/data/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneID")["geneSymbol"]
gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)["gene_name"]
c = c[c.sum(axis=1) >= 30]
c["sym"] = gi.reindex(c.index).values
c = c.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
nc = c / c.sum() * c.sum().median(); L = np.log2(nc + 1)
for line, pre in (("KP4", "KP4"), ("MiaPaCa2", "Mia")):
    a = L[[f"{pre}-T100-{i}" for i in (1, 2, 3)]]; b = L[[f"{pre}-DMSO-{i}" for i in (1, 2, 3)]]
    lf = a.mean(axis=1) - b.mean(axis=1); p = stats.ttest_ind(a, b, axis=1, equal_var=False).pvalue
    res[f"GSE206267 {line} trametinib 48h"] = pd.DataFrame({"lfc": lf, "p": p})
# ---- (2) GSE315130 FPKM
for line in ("PSN1", "SW1990"):
    d = pd.read_csv(f"pdac/GSE315130_{line}_cell_line_gene_sample_FPKM_with_symbol.txt.gz", sep="\t")
    d = d.drop_duplicates("Symbol").set_index("Symbol"); d = d[(d.filter(like="_Con_").mean(axis=1) > 0.5) | (d.filter(like="_6236_").mean(axis=1) > 0.5)]
    L2 = np.log2(d.filter(regex="_(Con|6236|Taze|Combo)_") + 1)
    for tag, lab in (("6236", "RMC-6236"), ("Combo", "RMC-6236 + tazemetostat")):
        a = L2[[f"{line}_{tag}_{i}" for i in (1, 2, 3)]]; b = L2[[f"{line}_Con_{i}" for i in (1, 2, 3)]]
        lf = a.mean(axis=1) - b.mean(axis=1); p = stats.ttest_ind(a, b, axis=1, equal_var=False).pvalue
        res[f"GSE315130 {line} {lab}"] = pd.DataFrame({"lfc": lf, "p": p})
# ---- add earlier datasets for the same table
for l in ("CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"):
    d = pd.read_csv(f"out/pdac_de_{l}_h24.tsv", sep="\t").set_index("gene")
    res[f"GSE238198 {l} trametinib 24h"] = pd.DataFrame({"lfc": d["logFC"], "p": d["P.Value"]})
for l, nm in (("PANC1", "PANC1_Simvastatin_High"), ("MiaPaca2", "MiaPaca2_Simvastatin_High")):
    d = pd.read_csv(f"out/statin_de_{nm}.tsv", sep="\t"); s = pd.read_csv("pdac/statin_symbols.tsv", sep="\t").set_index("gene_id")["gene_short_name"]
    d["sym"] = d.gene.map(s); d = d.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
    res[f"GSE149566 {l} simvastatin 20uM"] = pd.DataFrame({"lfc": d["logFC"], "p": d["P.Value"]})
print("[1] trunk genes (log2FC; * p<.05 unadjusted):")
hdr = f"{'dataset':42}" + "".join(f"{g:>9}" for g in TRUNK)
print(hdr)
for k, d in res.items():
    cells = []
    for g in TRUNK:
        if g in d.index: cells.append(f"{d.loc[g,'lfc']:+.2f}{'*' if d.loc[g,'p']<.05 else ' '}")
        else: cells.append("NA")
    print(f"{k:42}" + "".join(f"{x:>9}" for x in cells))
print("\n[2] program level, 24 listed SREBP2/sterol genes: n measured, n up (log2FC>0), mean log2FC | trunk (6 genes) mean | cholesterol branch (11 genes) mean")
rows = []
for k, d in res.items():
    s = d.reindex([g for g in SREBP if g in d.index])
    t = d.reindex([g for g in TRUNK if g in d.index]); b = d.reindex([g for g in BRANCH if g in d.index])
    bino = stats.binomtest(int((s.lfc > 0).sum()), len(s), 0.5).pvalue
    rows.append((k, len(s), int((s.lfc > 0).sum()), s.lfc.mean(), t.lfc.mean(), b.lfc.mean(), bino))
    print(f"  {k:42} n={len(s):2d} up={int((s.lfc>0).sum()):2d} mean {s.lfc.mean():+.2f} | trunk {t.lfc.mean():+.2f} | branch {b.lfc.mean():+.2f} | sign test p={bino:.3g}")
# background: how unusual is the program mean vs random 24-gene sets of expressed genes?
rng = np.random.default_rng(0)
print("\n[3] empirical check: percentile of the program mean among 2000 random 24-gene sets of the same dataset")
for k, d in res.items():
    s = d.reindex([g for g in SREBP if g in d.index]).lfc.dropna()
    pool = d.lfc.dropna().values
    rnd = np.array([rng.choice(pool, len(s), replace=False).mean() for _ in range(2000)])
    print(f"  {k:42} program mean {s.mean():+.2f}; random-set mean sd {rnd.std():.2f}; fraction of random sets >= observed {(rnd >= s.mean()).mean():.4f}")
