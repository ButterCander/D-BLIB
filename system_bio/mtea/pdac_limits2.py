"""Further limitation tests for the mevalonate-trunk (R2D) hypothesis:
(A) kinetics 3/12/24 h in four PDAC lines; (B) data-derived statin (SREBP2-like) signature outside the 24 listed genes;
(C) is the trunk induction just a correlate of cell-cycle suppression? (line x time points, plus within-line gene level)."""
import numpy as np, pandas as pd
from scipy import stats
rng = np.random.default_rng(1)
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TRUNK = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
S = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()
CC = S + G2M
LINES = ["CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"]
D = {(l, t): pd.read_csv(f"out/pdac_de_{l}_{t}.tsv", sep="\t").set_index("gene") for l in LINES for t in ("h3", "h12", "h24")}
print("[A] kinetics: mean log2FC of the 24 sterol-program genes | trunk (6) | cell-cycle signature (S+G2/M)")
rows = []
for l in LINES:
    out = []
    for t in ("h3", "h12", "h24"):
        d = D[(l, t)]["logFC"]
        p = d.reindex([g for g in SREBP if g in d.index]).mean(); tr = d.reindex([g for g in TRUNK if g in d.index]).mean(); cc = d.reindex([g for g in CC if g in d.index]).mean()
        out.append(f"{t}: {p:+.2f} | {tr:+.2f} | {cc:+.2f}"); rows.append((l, t, p, tr, cc))
    print(f"  {l:9} " + "   ".join(out))
R = pd.DataFrame(rows, columns=["line", "t", "prog", "trunk", "cc"])
rho, p = stats.spearmanr(R.prog, R.cc)
print(f"  across 12 line x time points: Spearman(program, cell-cycle) = {rho:+.2f} (p = {p:.3f}); Spearman(trunk, cell-cycle) = {stats.spearmanr(R.trunk, R.cc)[0]:+.2f}")
# (B) data-derived signature: genes induced by simvastatin in both PDAC lines (outside the listed 24)
st = {}
for l, nm in (("PANC1", "PANC1_Simvastatin_High"), ("MiaPaca2", "MiaPaca2_Simvastatin_High")):
    d = pd.read_csv(f"out/statin_de_{nm}.tsv", sep="\t"); s = pd.read_csv("pdac/statin_symbols.tsv", sep="\t").set_index("gene_id")["gene_short_name"]
    d["sym"] = d.gene.map(s); st[l] = d.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")
common = st["PANC1"].index.intersection(st["MiaPaca2"].index)
sc = pd.DataFrame({"a": st["PANC1"].loc[common, "logFC"], "b": st["MiaPaca2"].loc[common, "logFC"],
                   "pa": st["PANC1"].loc[common, "adj.P.Val"], "pb": st["MiaPaca2"].loc[common, "adj.P.Val"]})
sig = sc[(sc.a > 0.75) & (sc.b > 0.75) & (sc.pa < 0.05) & (sc.pb < 0.05)].copy()
sig["m"] = (sig.a + sig.b) / 2
sig = sig[~sig.index.isin(SREBP)].sort_values("m", ascending=False)
sigset = list(sig.index[:100])
print(f"\n[B] statin-derived signature outside the 24 listed genes: {len(sig)} genes up >0.75 in both lines; top-100 used. Examples: {', '.join(sigset[:15])}")
def score(d, genes):
    x = d.reindex(genes).dropna(); return x.mean(), len(x)
def pctl(d, genes, n=2000):
    x = d.reindex(genes).dropna(); pool = d.dropna().values
    r = np.array([rng.choice(pool, len(x), replace=False).mean() for _ in range(n)]); return (r >= x.mean()).mean()
print("  MEKi / statin response of this signature (mean log2FC; fraction of random 100-gene sets >= observed):")
for l in LINES:
    for t in ("h12", "h24"):
        d = D[(l, t)]["logFC"]; m, n = score(d, sigset); print(f"    GSE238198 {l:9} trametinib {t}: {m:+.2f} (n={n}; p~{pctl(d, sigset):.4f})")
# GSE206267 (48 h) recomputed here
c = pd.read_csv("pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz", sep="\t", index_col=0); c = c[c.sum(axis=1) >= 30]
gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)["gene_name"]; c["sym"] = gi.reindex(c.index).values
c = c.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym"); nc = c / c.sum() * c.sum().median(); L = np.log2(nc + 1)
for line, pre in (("KP4", "KP4"), ("MiaPaCa2", "Mia")):
    d = L[[f"{pre}-T100-{i}" for i in (1, 2, 3)]].mean(axis=1) - L[[f"{pre}-DMSO-{i}" for i in (1, 2, 3)]].mean(axis=1)
    m, n = score(d, sigset); print(f"    GSE206267 {line:9} trametinib 48h: {m:+.2f} (n={n}; p~{pctl(d, sigset):.4f})")
# AGS official Tsirvouli 24h
off = pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv"); gii = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)
for c_ in ("PD", "PI", "PIPD"):
    o = off[off.Contrast == f"Diff_{c_}_24h"].copy(); o["sym"] = gii.reindex(o.Gene_name)["gene_name"].values
    d = o.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")["logFC"]; m, n = score(d, sigset)
    print(f"    AGS Tsirvouli {c_:5} 24h: {m:+.2f} (n={n}; p~{pctl(d, sigset):.4f})")
# (C) within-line gene-level: is the induced set confounded with cell-cycle genes? overlap of sterol program with CC genes = 0 by construction; check correlation of per-sample?
print("\n[C] gene-level: among genes down-regulated with the cell cycle (CC signature), fraction of trunk genes that are also in CC list:", len(set(TRUNK) & set(CC)))
