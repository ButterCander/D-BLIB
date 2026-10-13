import json, os, time
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
t00 = time.time(); rng = np.random.default_rng(7)
OUT = "out"; D = "mteapy/src/mteapy/data"
# ---------------- data ----------------
cnt = pd.read_csv("GSE285616_counts.csv.gz", index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]; cnt = cnt.loc[cnt.sum(axis=1) >= 10]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
Y = np.log2(x / sf + 1)                       # genes x 24
G = Y.shape[0]; genes = list(cnt.index); gi = {g: i for i, g in enumerate(genes)}
batch = np.array([i // 6 for i in range(24)])
tasks = json.load(open(f"{OUT}/task_genes.json"))
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneID").set_index("geneID")["geneSymbol"]
ess = pd.read_csv(f"{D}/HumanGEM_essential_genes_matrix.tsv", sep="\t", index_col=0).astype(bool)
sets = {}
for t, gl in tasks.items():
    idx = sorted({gi[g] for g in gl if g in gi})
    if 5 <= len(idx) <= 200: sets["T:" + t] = np.array(idx)
for t in ess.columns:
    idx = sorted({gi[gm[g]] for g in ess.index[ess[t]] if g in gm.index and gm[g] in gi})
    if 3 <= len(idx) <= 200: sets["E:" + t] = np.array(idx)
names = list(sets); sizes = sorted({len(v) for v in sets.values()})
isT = np.array([n.startswith("T:") for n in names])
print(f"gene sets: {int(isT.sum())} task sets, {int((~isT).sum())} essential sets; genes {G}", flush=True)


def design(cond):
    return np.column_stack([np.ones(24)] + [(batch == b).astype(float) for b in (1, 2, 3)] + [(cond == c).astype(float) for c in range(1, 6)])


def fit(Yd, cond):
    X = design(cond); XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ X.T @ Yd.T
    R = Yd.T - X @ beta; s2 = (R ** 2).sum(0) / (24 - X.shape[1])
    Rz = (R - R.mean(0)).T; Rz /= (np.linalg.norm(Rz, axis=1, keepdims=True) + 1e-12)
    return beta, XtXi, s2, Rz


def zstat(beta, XtXi, s2, k):
    c = np.zeros(9); c[4 + k - 1] = 1
    est = c @ beta; se = np.sqrt(s2 * (c @ XtXi @ c)); t = est / se
    return stats.norm.ppf(np.clip(stats.t.cdf(t, 15), 1e-12, 1 - 1e-12))


def pvals(z, Rz, ndraw=500):
    n = len(z); o_ind = np.empty(len(names)); o_cam = np.empty(len(names)); o_perm = np.empty(len(names))
    null = {m: z[rng.integers(0, n, size=(ndraw, m))].mean(1) for m in sizes}
    for k, nm in enumerate(names):
        idx = sets[nm]; m = len(idx); inset = np.zeros(n, bool); inset[idx] = True
        zin = z[idx].mean(); zout = z[~inset]; so = zout.var(ddof=1)
        C = Rz[idx] @ Rz[idx].T; rho = max(0.0, (C.sum() - m) / (m * (m - 1)))
        for rr, store in ((0.0, o_ind), (rho, o_cam)):
            vif = 1 + (m - 1) * rr; t = (zin - zout.mean()) / np.sqrt(so * (vif / m + 1 / (n - m))); store[k] = 2 * stats.t.sf(abs(t), n - 2)
        nl = null[m]; o_perm[k] = (1 + np.sum(np.abs(nl - nl.mean()) >= abs(zin - nl.mean()))) / (ndraw + 1)
    return o_ind, o_cam, o_perm


# ---------------- A. consistency ----------------
print("\n[A] consistency of contrasts and conversions")
L = {}
for k, f in (("AB", "de_PI103+PD"), ("A", "de_PI103"), ("B", "de_PD0325901"), ("AB2", "de_PI103+5Z"), ("C", "de_5Z")):
    L[k] = pd.read_csv(f"{OUT}/{f}.tsv", sep="\t", index_col=0)["log2FoldChange"]


def chk(name, file, expr):
    r = pd.read_csv(f"{OUT}/{file}.tsv", sep="\t", index_col=0)["log2FoldChange"]; e = expr.reindex(r.index)
    ok = r.notna() & e.notna(); d = (r[ok] - e[ok]).abs()
    print(f"  {name:34} n={int(ok.sum())} max|diff|={d.max():.2e} corr={np.corrcoef(r[ok], e[ok])[0, 1]:.6f}")


chk("AB-A  (PI3Ki+MEKi vs PI3Ki)", "inter_PI3Ki_MEKi_vsPI3Ki", L["AB"] - L["A"])
chk("AB-B  (PI3Ki+MEKi vs MEKi)", "inter_PI3Ki_MEKi_vsMEKi", L["AB"] - L["B"])
chk("AB-.5A-.5B (dose-scaled)", "inter_PI3Ki_MEKi_dose", L["AB"] - .5 * L["A"] - .5 * L["B"])
chk("AB-A-B (additive)", "inter_PI3Ki_MEKi_add", L["AB"] - L["A"] - L["B"])
chk("AB2-C (PI3Ki+TAKi vs TAKi)", "inter_PI3Ki_TAKi_vsTAKi", L["AB2"] - L["C"])
a = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki-MEKi", header=1).iloc[:, :5]; a.columns = ["t", "d", "s", "p", "sig"]; a = a.dropna(subset=["t"]); a["p"] = pd.to_numeric(a.p)
p2 = np.minimum(1, 2 * a.p.values); q = false_discovery_control(p2, method="bh")


def bh_manual(p):
    o = np.argsort(p); n = len(p); r = p[o] * n / np.arange(1, n + 1); r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(r, 1); return out


qt = false_discovery_control(a.p.values, method="bh")
print(f"  BH scipy vs manual: max|diff|={np.abs(q - bh_manual(p2)).max():.2e}; authors PI3Ki-MEKi q<.05 = {int((q < .05).sum())}")
print(f"  p convention (authors, PI3Ki-MEKi, 187 tasks): two-sided 2p & q<.05: {int((q < .05).sum())} | tail p used directly & q<.05: {int((qt < .05).sum())} | tail p & q<.025: {int((qt < .025).sum())}")

# ---------------- B. null calibration ----------------
print("\n[B] negative control: condition labels shuffled within batch (all drug effects removed)", flush=True)
NS = 30
allp = {f: {k: [] for k in ("ind", "cam", "perm")} for f in "TE"}; nrej = {f: {k: [] for k in ("ind", "cam", "perm")} for f in "TE"}
for s in range(NS):
    cond = np.concatenate([rng.permutation(6) for _ in range(4)])
    beta, XtXi, s2, Rz = fit(Y, cond)
    for k in range(1, 6):
        z = zstat(beta, XtXi, s2, k); pi, pc, pp = pvals(z, Rz)
        for fam, mask in (("T", isT), ("E", ~isT)):
            for key, pv in (("ind", pi), ("cam", pc), ("perm", pp)):
                v = pv[mask]; allp[fam][key].append(v); nrej[fam][key].append(int((false_discovery_control(v, method="bh") < .05).sum()))
    if s % 10 == 9: print(f"   shuffles done: {s + 1}/{NS} ({time.time() - t00:.0f}s)", flush=True)
print("  false-positive rate and mean BH rejections per null contrast (150 contrasts; nominal 0.05 / 0):")
for fam, lab in (("T", f"task gene sets (n={int(isT.sum())})"), ("E", f"essential-gene sets (n={int((~isT).sum())})")):
    for key, lab2 in (("ind", "parametric, independence"), ("cam", "CAMERA-like, estimated rho"), ("perm", "gene-label permutation (as MTEApy)")):
        v = np.concatenate(allp[fam][key]); ks = stats.kstest(v, "uniform").statistic
        print(f"   {lab:28} {lab2:36} FPR@.05={np.mean(v < .05):.3f} FPR@.01={np.mean(v < .01):.3f} KS={ks:.3f} meanBH={np.mean(nrej[fam][key]):.2f}", flush=True)

# ---------------- C. positive control ----------------
print("\n[C] positive control: add a known log2 shift to one task's genes in the 4 samples of one condition", flush=True)
tnames = [n for n in names if n.startswith("T:")]; tpos = {n: i for i, n in enumerate(tnames)}
for delta in (0.1, 0.2, 0.3, 0.5):
    det = {"ind": 0, "cam": 0, "perm": 0}; bh = {"ind": 0, "cam": 0, "perm": 0}; trials = 60
    for _ in range(trials):
        cond = np.concatenate([rng.permutation(6) for _ in range(4)])
        nm = tnames[rng.integers(len(tnames))]; Ys = Y.copy(); Ys[np.ix_(sets[nm], np.where(cond == 1)[0])] += delta
        beta, XtXi, s2, Rz = fit(Ys, cond); z = zstat(beta, XtXi, s2, 1); pi, pc, pp = pvals(z, Rz)
        k = names.index(nm)
        for key, pv in (("ind", pi), ("cam", pc), ("perm", pp)):
            det[key] += pv[k] < .05
            qq = false_discovery_control(pv[isT], method="bh"); bh[key] += qq[tpos[nm]] < .05
    print(f"   delta={delta:.1f} log2: p<.05 detect: indep {det['ind'] / trials:.2f} | CAMERA-like {det['cam'] / trials:.2f} | permutation {det['perm'] / trials:.2f} ; BH q<.05: indep {bh['ind'] / trials:.2f} | CAMERA-like {bh['cam'] / trials:.2f} | permutation {bh['perm'] / trials:.2f}", flush=True)

# ---------------- D. detectable difference for the HSA contrast ----------------
print("\n[D] smallest task-score magnitude the HSA contrast (AB vs PI3Ki) could call, by set size")
r = pd.read_csv(f"{OUT}/inter_PI3Ki_MEKi_vsPI3Ki.tsv", sep="\t", index_col=0)
v = np.where(r.padj.fillna(1) < .05, r.log2FoldChange.fillna(0), 0.0)
for m in (6, 17, 50):
    mm = v[rng.integers(0, len(v), size=(20000, m))].mean(1)
    print(f"   set size {m:3d}: null SD of mean = {mm.std():.3f}; |score| needed for p<.05 two-sided = {np.quantile(np.abs(mm), .95):.3f}; for p<.001 = {np.quantile(np.abs(mm), .999):.3f}")
print("   authors weak-synergy margins (median) = 0.076-0.081 log2FC")
print(f"\ntotal {time.time() - t00:.0f}s")
