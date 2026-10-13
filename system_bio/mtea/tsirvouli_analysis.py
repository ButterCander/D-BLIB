"""Independent, time-resolved RNA-seq of the same AGS system (Tsirvouli et al. 2024 repo): log2-normalised expression,
6 treatments x 5 time points x 3 replicates. Checks: dosing design, replication of the 24 h result, HSA / interaction at gene-set level, proliferation."""
import json, re
import numpy as np, pandas as pd
from scipy import stats, special
from scipy.stats import false_discovery_control
rng = np.random.default_rng(3)
Y = pd.read_csv("tsirvouli/counts_norm.csv", index_col=0)
gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)
sym = gi.reindex(Y.index)["gene_name"]
ok = sym.notna() & ~sym.duplicated()
Y = Y[ok.values]; Y.index = sym[ok].values
meta = pd.DataFrame([re.match(r"([A-Za-z0-9]+)_(\d+)_Rep(\d)", c).groups() for c in Y.columns], columns=["trt", "t", "rep"], index=Y.columns)
meta["t"] = meta.t.astype(int)
TRT = ["DMSO", "OXO", "PD", "PI", "PIOXO", "PIPD"]; TIMES = [1, 2, 4, 8, 24]
print(f"genes {Y.shape[0]}, samples {Y.shape[1]}", flush=True)

# ---------------- limma-like moderated t ----------------
def trigamma_inv(x):
    y = 0.5 + 1 / x
    for _ in range(50):
        y = y + (special.polygamma(1, y) * (1 - special.polygamma(1, y) / x)) / special.polygamma(2, y)
    return y
def fit_t(Ysub, trt, contrasts):
    """cell-means model; returns dict contrast -> (est, moderated z) for every gene"""
    lev = TRT; n = Ysub.shape[1]
    M = np.array([[1.0 if trt[j] == l else 0.0 for l in lev] for j in range(n)])
    B = np.linalg.lstsq(M, Ysub.T, rcond=None)[0]           # 6 x G
    R = Ysub.T - M @ B; d = n - len(lev); s2 = (R ** 2).sum(0) / d
    z = np.log(s2 + 1e-12); e = z - special.digamma(d / 2) + np.log(d / 2)
    em = e.mean(); ev = e.var(ddof=1) - special.polygamma(1, d / 2)
    if ev > 0:
        d0 = 2 * trigamma_inv(ev); s0 = np.exp(em + special.digamma(d0 / 2) - np.log(d0 / 2))
    else:
        d0 = np.inf; s0 = np.exp(em)
    s2t = (d0 * s0 + d * s2) / (d0 + d) if np.isfinite(d0) else np.full_like(s2, s0)
    dft = d + (d0 if np.isfinite(d0) else 1e6)
    cnt = M.sum(0); out = {}
    for name, w in contrasts.items():
        c = np.array([w.get(l, 0.0) for l in lev]); est = c @ B; v = (c ** 2 / cnt).sum()
        t = est / np.sqrt(s2t * v)
        out[name] = (est, stats.norm.ppf(np.clip(stats.t.cdf(t, dft), 1e-12, 1 - 1e-12)), 2 * stats.t.sf(np.abs(t), dft))
    return out, d0
CON = {
    "PI": {"PI": 1, "DMSO": -1}, "PD": {"PD": 1, "DMSO": -1}, "OXO": {"OXO": 1, "DMSO": -1},
    "PIPD": {"PIPD": 1, "DMSO": -1}, "PIOXO": {"PIOXO": 1, "DMSO": -1},
    "PIPD-PI": {"PIPD": 1, "PI": -1}, "PIPD-PD": {"PIPD": 1, "PD": -1},
    "PIOXO-PI": {"PIOXO": 1, "PI": -1}, "PIOXO-OXO": {"PIOXO": 1, "OXO": -1},
    "PIPD-int": {"PIPD": 1, "PI": -1, "PD": -1, "DMSO": 1}, "PIOXO-int": {"PIOXO": 1, "PI": -1, "OXO": -1, "DMSO": 1},
}
RES = {}
for t in TIMES:
    cols = meta.index[(meta.t == t) & meta.trt.isin(TRT)]
    sub = Y[cols].values; trt = meta.loc[cols, "trt"].tolist()
    RES[t], d0 = fit_t(sub, trt, CON)
    print(f"t={t:2d}h: {len(cols)} samples, prior df d0={d0:.1f}", flush=True)
genes = list(Y.index)

# ---------------- 1. dosing design diagnostic ----------------
print("\n[1] Is the combination dosed like the singles? slope of LFC(combination vs DMSO) on LFC(single vs DMSO)")
print("    genes with |LFC single|>0.5 and |LFC partner|<0.25; slope through origin; 1.0 = same dose, ~0.5-0.7 = reduced dose")
for t in TIMES:
    row = []
    for combo, a, b in (("PIPD", "PI", "PD"), ("PIOXO", "PI", "OXO")):
        for single, other in ((a, b), (b, a)):
            la, lo, lc = RES[t][single][0], RES[t][other][0], RES[t][combo][0]
            m = (np.abs(la) > 0.5) & (np.abs(lo) < 0.25)
            if m.sum() < 20: row.append(f"{combo}~{single}: n/a"); continue
            row.append(f"{combo}~{single}: {np.sum(la[m]*lc[m])/np.sum(la[m]**2):.2f} (n={int(m.sum())})")
    print(f"    {t:2d}h  " + " | ".join(row))
de = {k: pd.read_csv(f"out/de_{k}.tsv", sep="\t", index_col=0)["log2FoldChange"] for k in ("PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD")}
for combo, a, b in (("PI103+PD", "PI103", "PD0325901"), ("PI103+5Z", "PI103", "5Z")):
    for single, other in ((a, b), (b, a)):
        la = de[single]; lo = de[other].reindex(la.index); lc = de[combo].reindex(la.index)
        m = (la.abs() > 0.5) & (lo.abs() < 0.25)
        print(f"    GSE285616 24h  {combo}~{single}: {np.sum(la[m]*lc[m])/np.sum(la[m]**2):.2f} (n={int(m.sum())})")

# ---------------- 2. replication of the 24 h gene-level result ----------------
print("\n[2] replication at 24 h: log2FC vs DMSO in the independent experiment vs GSE285616 (per-gene correlation, shared genes)")
mapn = {"PI": "PI103", "PD": "PD0325901", "OXO": "5Z", "PIPD": "PI103+PD", "PIOXO": "PI103+5Z"}
g_index = pd.Index(genes)
for k, g in mapn.items():
    ref = de[g]; common = ref.index.intersection(g_index); ref = ref.loc[common]
    mine = pd.Series(RES[24][k][0], index=g_index).loc[common]
    qv = false_discovery_control(RES[24][k][2], method="bh"); ndeg = int(((qv < .05) & (np.abs(RES[24][k][0]) >= 1)).sum())
    print(f"    {k:6} vs {g:10}: n={len(common)}  Pearson r={np.corrcoef(ref, mine)[0,1]:.3f}  Spearman={stats.spearmanr(ref, mine)[0]:.3f} | DEG(padj<.05,|lfc|>=1) here {ndeg}")
key = ["ASNS", "PHGDH", "PSAT1", "PSPH", "SHMT2", "ODC1", "SRM", "AMD1", "TXNRD1", "COQ3", "UPB1"]
print("    key genes log2FC vs DMSO (24 h): " + " | ".join(f"{k}: " + ", ".join(f"{c} {RES[24][c][0][genes.index(k)]:+.1f}" for c in ("PI", "PD", "OXO", "PIPD", "PIOXO")) for k in key[:4] if k in genes))

# ---------------- 3. gene-set (task) level HSA / interaction, time resolved ----------------
tasks = json.load(open("out/task_genes.json"))
gdict = {g: i for i, g in enumerate(genes)}
sets = {t: np.array(sorted({gdict[g] for g in gl if g in gdict})) for t, gl in tasks.items()}
sets = {t: v for t, v in sets.items() if 5 <= len(v) <= 200}
# residual correlation from the full cell-means model over all 94 samples
cell = meta.trt + "_" + meta.t.astype(str)
Yc = Y.values - pd.DataFrame(Y.values.T, index=Y.columns).groupby(cell.values).transform("mean").values.T
Rz = Yc - Yc.mean(1, keepdims=True); Rz /= (np.linalg.norm(Rz, axis=1, keepdims=True) + 1e-12)
rho = {t: max(0.0, (lambda C, m: (C.sum() - m) / (m * (m - 1)))(Rz[v] @ Rz[v].T, len(v))) for t, v in sets.items()}
print(f"\n[3] gene-set (task) tests, {len(sets)} task sets, CAMERA-like with estimated rho (median rho {np.median(list(rho.values())):.2f})")
def settest(z):
    n = len(z); rows = []
    for t, idx in sets.items():
        inset = np.zeros(n, bool); inset[idx] = True; m = len(idx)
        zin = z[idx].mean(); zout = z[~inset]; vif = 1 + (m - 1) * rho[t]
        tt = (zin - zout.mean()) / np.sqrt(zout.var(ddof=1) * (vif / m + 1 / (n - m)))
        rows.append((t, zin, 2 * stats.t.sf(abs(tt), n - 2)))
    d = pd.DataFrame(rows, columns=["task", "z", "p"]).set_index("task"); d["q"] = false_discovery_control(d.p.values, method="bh")
    return d
hdr = f"{'contrast':11}" + "".join(f"{str(t)+'h':>16}" for t in TIMES)
print("    entries: p<.05 (neg) [BH<.05]   -- negative = more down-regulated")
print("    " + hdr)
for name in ("PI", "PD", "OXO", "PIPD", "PIOXO", "PIPD-PI", "PIPD-PD", "PIOXO-PI", "PIOXO-OXO", "PIPD-int", "PIOXO-int"):
    cells = []
    for t in TIMES:
        d = settest(RES[t][name][1]); s = d[d.p < .05]
        cells.append(f"{len(s):3d} ({int((s.z<0).sum()):3d}) [{int((d.q<.05).sum()):3d}]")
    print(f"    {name:11}" + "".join(f"{c:>16}" for c in cells))

# ---------------- 4. proliferation signature over time ----------------
S = "MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 CENPU HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8".split()
G2M = "HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 PIMREG SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 JPT1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA".split()
idx = [gdict[g] for g in S + G2M if g in gdict]
print(f"\n[4] cell-cycle signature ({len(idx)} genes), mean log2FC vs DMSO at the same time point")
print("    " + f"{'':8}" + "".join(f"{str(t)+'h':>9}" for t in TIMES))
for name in ("PI", "PD", "OXO", "PIPD", "PIOXO"):
    print(f"    {name:8}" + "".join(f"{RES[t][name][0][idx].mean():>+9.2f}" for t in TIMES))
print("    growth at 48 h (Flobak 2015, digitized): PI 52-73% | PD 47% | 5Z 59% | PI+PD 78% | PI+5Z 92% inhibited")

# ---------------- 5. which tasks carry the few significant signals ----------------
meta_t = pd.read_csv("mteapy/src/mteapy/data/task_metadata.tsv", sep="\t").set_index("ID")
print("\n[5] top task sets for selected contrasts (z<0 = more down-regulated; q = BH over task sets)")
for name, t in (("PIPD-int", 8), ("PIPD-int", 24), ("PIPD-PI", 24), ("PIOXO-PI", 8), ("PIOXO-int", 24), ("PIPD", 24), ("PI", 24)):
    d = settest(RES[t][name][1]).sort_values("p").head(4)
    print(f"    {name:9} {t:2d}h: " + " | ".join(f"{i} {meta_t.loc[i,'DESCRIPTION'][:26]} z={r.z:+.2f} q={r.q:.3f}" for i, r in d.iterrows()))
