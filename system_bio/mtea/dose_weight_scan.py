"""Which half-dose weight w (expected combination effect = w*A + w*B) makes the metabolic 'synergy' disappear?
Growth dose-response curves (Flobak 2015 S1) put the half-GI50 effect at about 0.7-0.8 of the full-GI50 effect."""
import json
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
OUT = "out"; D = "mteapy/src/mteapy/data"
cnt = pd.read_csv("GSE285616_counts.csv.gz", index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]; cnt = cnt.loc[cnt.sum(axis=1) >= 10]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
Y = np.log2(x / sf + 1)
batch = np.array([i // 6 for i in range(24)]); cond = np.array([i % 6 for i in range(24)])
X = np.column_stack([np.ones(24)] + [(batch == b).astype(float) for b in (1, 2, 3)] + [(cond == c).astype(float) for c in range(1, 6)])
XtXi = np.linalg.inv(X.T @ X); beta = XtXi @ X.T @ Y.T
R = Y.T - X @ beta; s2 = (R ** 2).sum(0) / 15
Rz = (R - R.mean(0)).T; Rz /= (np.linalg.norm(Rz, axis=1, keepdims=True) + 1e-12)
genes = list(cnt.index); gi = {g: i for i, g in enumerate(genes)}
tasks = json.load(open(f"{OUT}/task_genes.json"))
sets = {t: np.array(sorted({gi[g] for g in gl if g in gi})) for t, gl in tasks.items()}
sets = {t: v for t, v in sets.items() if 5 <= len(v) <= 200}
rho = {}
for t, v in sets.items():
    C = Rz[v] @ Rz[v].T; m = len(v); rho[t] = max(0.0, (C.sum() - m) / (m * (m - 1)))
syn = pd.ExcelFile("si/MOESM5.xlsx"); auth = {}
for combo in ("PI3Ki-TAKi", "PI3Ki-MEKi"):
    d = syn.parse(combo).iloc[:, :2]; auth[combo] = {("X" + str(i).lstrip("X")) for i in d.iloc[:, 0].dropna().astype(str)}
COEF = {"PI103": 4, "PD0325901": 5, "5Z": 6, "PI103+5Z": 7, "PI103+PD": 8}


def run(combo_key, wa, wb):
    a, b, ab = (("PI103", "PD0325901", "PI103+PD") if combo_key == "MEKi" else ("PI103", "5Z", "PI103+5Z"))
    c = np.zeros(9); c[COEF[ab]] = 1; c[COEF[a]] = -wa; c[COEF[b]] = -wb
    est = c @ beta; se = np.sqrt(s2 * (c @ XtXi @ c)); t = est / se
    z = stats.norm.ppf(np.clip(stats.t.cdf(t, 15), 1e-12, 1 - 1e-12)); n = len(z)
    rows = []
    for tk, idx in sets.items():
        inset = np.zeros(n, bool); inset[idx] = True; m = len(idx)
        zin = z[idx].mean(); zout = z[~inset]; vif = 1 + (m - 1) * rho[tk]
        tt = (zin - zout.mean()) / np.sqrt(zout.var(ddof=1) * (vif / m + 1 / (n - m)))
        rows.append((tk, zin, 2 * stats.t.sf(abs(tt), n - 2)))
    d = pd.DataFrame(rows, columns=["task", "z", "p"]).set_index("task"); d["q"] = false_discovery_control(d.p.values, method="bh")
    return d


print("Interaction test vs assumed half-dose weight w  (expected AB = wA*A + wB*B); CAMERA-like, 120 task sets")
print(f"{'combo':12}{'wA/wB':>12}{'p<.05':>8}{'neg':>6}{'BH<.05':>8}{'BH<.10':>8}{'authors calls hit (p<.05)':>28}")
for ck, name in (("MEKi", "PI3Ki-MEKi"), ("TAKi", "PI3Ki-TAKi")):
    for wa, wb in ((0.5, 0.5), (0.6, 0.6), (0.7, 0.7), (0.7, 0.8), (0.8, 0.8), (0.9, 0.9), (1.0, 1.0)):
        d = run(ck, wa, wb)
        sig = d[d.p < .05]
        print(f"{name:12}{wa:>6.1f}/{wb:<5.1f}{len(sig):>8}{int((sig.z < 0).sum()):>6}{int((d.q < .05).sum()):>8}{int((d.q < .10).sum()):>8}{len(set(sig.index) & auth[name]):>14} of {len(auth[name] & set(d.index))}")
