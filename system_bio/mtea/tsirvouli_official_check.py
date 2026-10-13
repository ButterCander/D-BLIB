"""Re-check the Tsirvouli-based statements with the authors' official limma-voom logFC (DEGs_perContrast_all.csv)
instead of VSN-scale differences. Compares the two scales and recomputes: dose slope, 24 h replication, key genes,
cell-cycle signature, nucleotide / serine / polyamine gene groups."""
import re
import numpy as np, pandas as pd
from scipy import stats
off = pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv")
gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)
off["sym"] = gi.reindex(off.Gene_name)["gene_name"].values
off = off.dropna(subset=["sym"]).drop_duplicates(["sym", "Contrast"])
LF = off.pivot(index="sym", columns="Contrast", values="logFC")
Y = pd.read_csv("tsirvouli/counts_norm.csv", index_col=0)
sym = gi.reindex(Y.index)["gene_name"]; ok = sym.notna() & ~sym.duplicated(); Y = Y[ok.values]; Y.index = sym[ok].values
meta = pd.DataFrame([re.match(r"([A-Za-z0-9]+)_(\d+)_Rep(\d)", c).groups() for c in Y.columns], columns=["trt", "t", "rep"], index=Y.columns); meta["t"] = meta.t.astype(int)
TREATS = ["PI", "PD", "OXO", "PIPD", "PIOXO"]; TIMES = [1, 2, 4, 8, 24]
def vsn(trt, t):
    a = Y[meta.index[(meta.trt == trt) & (meta.t == t)]].mean(axis=1); b = Y[meta.index[(meta.trt == "DMSO") & (meta.t == t)]].mean(axis=1)
    return a - b
V = {(c, t): vsn(c, t) for c in TREATS for t in TIMES}
O = {(c, t): LF[f"Diff_{c}_{t}h"].dropna() for c in TREATS for t in TIMES}
common = LF.index.intersection(Y.index)
print(f"official genes {LF.shape[0]}; VSN genes {Y.shape[0]}; shared {len(common)}")

print("\n[0] scale comparison, official logFC vs VSN-scale difference (shared genes)")
mean_expr = Y.mean(axis=1).reindex(common)
low = mean_expr < 7
rows = []
for c, t in ((c, t) for c in TREATS for t in (8, 24)):
    o = O[(c, t)].reindex(common); v = V[(c, t)].reindex(common); m = o.notna() & v.notna()
    slope = np.sum(o[m] * v[m]) / np.sum(o[m] ** 2)
    sl_low = np.sum(o[m & low] * v[m & low]) / np.sum(o[m & low] ** 2); sl_hi = np.sum(o[m & ~low] * v[m & ~low]) / np.sum(o[m & ~low] ** 2)
    rows.append(f"{c:5} {t:2d}h r={np.corrcoef(o[m], v[m])[0,1]:.3f} slope VSN~official all {slope:.2f} | mean<7 {sl_low:.2f} | mean>=7 {sl_hi:.2f}")
print("    " + "\n    ".join(rows))

print("\n[1] dose slope (genes with |single|>0.5, |partner|<0.25): official logFC | VSN diff")
for t in TIMES:
    out = []
    for combo, a, b in (("PIPD", "PI", "PD"), ("PIOXO", "PI", "OXO")):
        for single, other in ((a, b), (b, a)):
            res = []
            for S in (O, V):
                la, lo, lc = S[(single, t)], S[(other, t)], S[(combo, t)]
                g = la.index.intersection(lo.index).intersection(lc.index); la, lo, lc = la[g], lo[g], lc[g]
                m = (la.abs() > 0.5) & (lo.abs() < 0.25)
                res.append(f"{np.sum(la[m]*lc[m])/np.sum(la[m]**2):.2f}(n={int(m.sum())})" if m.sum() >= 20 else "n/a")
            out.append(f"{combo}~{single}: {res[0]} | {res[1]}")
    print(f"    {t:2d}h  " + "  ;  ".join(out))

print("\n[2] 24 h replication vs GSE285616 log2FC (shared genes), official | VSN")
de = {k: pd.read_csv(f"out/de_{k}.tsv", sep="\t", index_col=0)["log2FoldChange"] for k in ("PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD")}
mapn = {"PI": "PI103", "PD": "PD0325901", "OXO": "5Z", "PIPD": "PI103+PD", "PIOXO": "PI103+5Z"}
for k, g in mapn.items():
    ref = de[g].dropna(); res = []
    for S in (O, V):
        c = ref.index.intersection(S[(k, 24)].index); res.append(f"n={len(c)} r={np.corrcoef(ref[c], S[(k, 24)][c])[0,1]:.2f} rho={stats.spearmanr(ref[c], S[(k, 24)][c])[0]:.2f}")
    print(f"    {k:6} vs {g:10}: official {res[0]} | VSN {res[1]}")

print("\n[3] key genes, 24 h logFC (PI, PD, OXO, PIPD, PIOXO): official | VSN")
for g in ["ASNS", "PHGDH", "PSAT1", "PSPH", "SHMT2", "ODC1", "SRM", "AMD1", "TXNRD1", "COQ3"]:
    if g not in LF.index: continue
    print(f"    {g:6} official " + " ".join(f"{O[(c,24)].get(g, np.nan):+.2f}" for c in TREATS) + " | VSN " + " ".join(f"{V[(c,24)].get(g, np.nan):+.2f}" for c in TREATS))

print("\n[4] cell-cycle signature (S + G2/M), mean logFC vs DMSO: official | VSN")
S_ = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()
sig = [g for g in S_ + G2M if g in common]
print(f"    {len(sig)} genes;   " + "".join(f"{str(t)+'h':>14}" for t in TIMES))
for c in TREATS:
    print(f"    {c:6}" + "".join(f"{O[(c,t)].reindex(sig).mean():>+7.2f}|{V[(c,t)].reindex(sig).mean():+.2f} " for t in TIMES))

print("\n[5] metabolic gene-group mean logFC (8 h and 24 h): official | VSN")
G = {"purine de novo": "PPAT GART PFAS PAICS ADSL ATIC ADSS2 IMPDH2 GMPS".split(), "pyrimidine": "CAD DHODH UMPS CTPS1".split(),
     "dNTP/RNR-TYMS": "RRM1 RRM2 TYMS DUT TK1 DHFR".split(), "serine-glycine": "PHGDH PSAT1 PSPH SHMT2 MTHFD2".split(), "polyamine": "ODC1 SRM SMS AMD1".split()}
for name, gl in G.items():
    gl = [g for g in gl if g in common]
    print(f"    {name:16} n={len(gl)}")
    for t in (8, 24):
        print(f"        {t:2d}h " + "  ".join(f"{c} {O[(c,t)].reindex(gl).mean():+.2f}|{V[(c,t)].reindex(gl).mean():+.2f}" for c in TREATS))
