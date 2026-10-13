"""CAMERA-like competitive gene-set test (Wu & Smyth 2012) on the task gene sets.
Statistic: mean DESeq2 Wald z of the set minus mean z of the other genes; variance inflated by VIF = 1 + (m-1) * rho,
rho = fixed 0.01 (limma default) or the set's average pairwise correlation of expression residuals (clipped at 0)."""
import json, os
import numpy as np, pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = os.path.join(HERE, "out"); D = os.path.join(HERE, "mteapy", "src", "mteapy", "data")
CONDS = ["control", "PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]; TREATS = CONDS[1:]
cnt = pd.read_csv(os.path.join(HERE, "GSE285616_counts.csv.gz"), index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]
cnt = cnt.loc[cnt.sum(axis=1) >= 10]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
L = np.log2(x / sf + 1)
cond = np.array([i % 6 for i in range(24)]); batch = np.array([i // 6 for i in range(24)])
X = np.column_stack([np.ones(24)] + [(batch == b).astype(float) for b in (1, 2, 3)] + [(cond == c).astype(float) for c in range(1, 6)])
beta, *_ = np.linalg.lstsq(X, L.T, rcond=None); R = (L.T - X @ beta).T          # genes x samples residuals
Rz = (R - R.mean(1, keepdims=True)); Rz /= (np.linalg.norm(Rz, axis=1, keepdims=True) + 1e-12)
sym = {g: i for i, g in enumerate(cnt.index)}
tasks = json.load(open(os.path.join(OUT, "task_genes.json")))
gm = pd.read_csv(os.path.join(D, "HumanGEM_genes.tsv"), sep="\t").drop_duplicates("geneID").set_index("geneID")["geneSymbol"]
ess = pd.read_csv(os.path.join(D, "HumanGEM_essential_genes_matrix.tsv"), sep="\t", index_col=0).astype(bool)
ess_sets = {t: [gm[g] for g in ess.index[ess[t]] if g in gm.index] for t in ess.columns}

def rho_hat(genes):
    idx = [sym[g] for g in genes if g in sym]
    if len(idx) < 2: return 0.0
    C = Rz[idx] @ Rz[idx].T; m = len(idx)
    return max(0.0, (C.sum() - m) / (m * (m - 1)))

def camera(z, genes, rho):
    idx = np.array([zi for zi, g in zip(range(len(z)), z.index) if g in set(genes)])
    inset = z.index.isin(genes); m = int(inset.sum()); n = len(z)
    if m < 2 or m >= n - 1: return np.nan, np.nan, np.nan
    zv = z.values; zin, zout = zv[inset], zv[~inset]
    vif = 1 + (m - 1) * rho
    s2 = zout.var(ddof=1)
    t = (zin.mean() - zout.mean()) / np.sqrt(s2 * (vif / m + 1 / (n - m)))
    return t, 2 * stats.t.sf(abs(t), n - 2), vif

rows = []
for t in TREATS:
    d = pd.read_csv(os.path.join(OUT, f"de_{t}.tsv"), sep="\t", index_col=0).dropna(subset=["stat"])
    z = d["stat"]
    for fam, sets in (("TIDE-e genes", ess_sets), ("task genes", tasks)):
        for tid, genes in sets.items():
            gg = [g for g in genes if g in z.index]
            if len(gg) < 2: continue
            r_hat = rho_hat(gg)
            t0, p0, v0 = camera(z, gg, 0.0)
            t1, p1, v1 = camera(z, gg, 0.01)
            t2, p2, v2 = camera(z, gg, r_hat)
            rows.append(dict(cond=t, family=fam, task=tid, m=len(gg), rho=r_hat, t_indep=t0, p_indep=p0, p_fixed=p1, p_est=p2, vif_est=v2))
res = pd.DataFrame(rows); res.to_csv(os.path.join(OUT, "camera_like.tsv"), sep="\t", index=False)
print(res.groupby("family").agg(n=("task", "size"), median_m=("m", "median"), median_rho=("rho", "median"), q90_rho=("rho", lambda s: s.quantile(.9)), median_vif=("vif_est", "median")))
