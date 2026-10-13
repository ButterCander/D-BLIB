"""Does the 24 h RNA-seq carry the growth effect? Cell-cycle (Tirosh S + G2/M) signature per condition and for AB vs best single."""
import numpy as np, pandas as pd
from scipy import stats
S = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()
cnt = pd.read_csv("GSE285616_counts.csv.gz", index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]; cnt = cnt.loc[cnt.sum(axis=1) >= 10]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
Y = np.log2(x / sf + 1); genes = list(cnt.index)
cond = np.array([i % 6 for i in range(24)]); batch = np.array([i // 6 for i in range(24)])
names = ["control", "PI3Ki", "MEKi", "TAKi", "PI3Ki+TAKi", "PI3Ki+MEKi"]
X = np.column_stack([np.ones(24)] + [(batch == b).astype(float) for b in (1, 2, 3)] + [(cond == c).astype(float) for c in range(1, 6)])
XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ X.T @ Y.T; R = Y.T - X @ beta; s2 = (R ** 2).sum(0) / 15
for label, gl in (("S-phase", S), ("G2/M", G2M), ("S + G2/M", S + G2M)):
    idx = [genes.index(g) for g in gl if g in genes]
    print(f"\n{label}: {len(idx)} genes present")
    # condition effects vs control (mean log2FC of the signature genes), batch-adjusted
    eff = {names[k]: beta[3 + k, idx].mean() for k in range(1, 6)}
    print("  mean log2FC vs control:", {k: round(v, 3) for k, v in eff.items()})
    # contrasts with set-level t (genes treated as independent -> optimistic; plus replicate-level paired view)
    def contrast(cvec):
        c = np.zeros(9)
        for i, w in cvec.items(): c[i] = w
        est = (c @ beta)[idx]; se = np.sqrt(s2[idx] * (c @ XtXi @ c)); t = est / se
        return est.mean(), t.mean()
    for lab, cv in (("PI3Ki+MEKi vs PI3Ki alone", {8: 1, 4: -1}), ("PI3Ki+MEKi vs MEKi alone", {8: 1, 5: -1}),
                    ("PI3Ki+TAKi vs PI3Ki alone", {7: 1, 4: -1}), ("PI3Ki+TAKi vs TAKi alone", {7: 1, 6: -1}),
                    ("PI3Ki vs MEKi", {4: 1, 5: -1})):
        e, t = contrast(cv)
        print(f"  {lab:28} mean log2FC {e:+.3f}   mean t {t:+.2f}")
    # replicate-level: signature score per sample (mean of gene-wise centred log2 expression), batch-corrected by subtracting batch mean
    M = Y[idx]; M = M - M.mean(1, keepdims=True)
    sc = M.mean(0)
    bm = np.array([sc[batch == b].mean() for b in range(4)])
    adj = sc - bm[batch]
    print("  replicate signature (batch-centred), per condition mean ± sd:", {names[k]: f"{adj[cond == k].mean():+.2f}±{adj[cond == k].std(ddof=1):.2f}" for k in range(6)})
print("\nGrowth at 48 h after drug (digitized from Flobak 2015 S1 Fig 18 and 21; increment over control, drug added at 1.0):")
print("  PI3Ki+MEKi experiment: control +3.25 | PD 35 nM +1.71 (47% inhibited) | PI103 0.7 uM +1.55 (52%) | combination +0.73 (78%)")
print("  PI3Ki+TAKi experiment: control +3.50 | PI103 +0.96 (73%) | 5Z 0.5 uM +1.45 (59%) | combination +0.27 (92%)")
