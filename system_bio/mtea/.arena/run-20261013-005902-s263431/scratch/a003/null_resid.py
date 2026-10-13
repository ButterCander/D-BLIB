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



cond_true=np.array([i%6 for i in range(24)])
Xt=design(cond_true); beta_t=np.linalg.lstsq(Xt,Y.T,rcond=None)[0]
R=(Y.T-Xt@beta_t).T   # genes x 24 residuals: batch + true condition effects removed
Yres=R+Y.mean(axis=1,keepdims=True)
NS=30
res={}
for lab,YY in (("original data (drug effects remain in residuals)",Y),("residual data (batch+condition removed = true null)",Yres)):
    allp={f:{k:[] for k in ("ind","cam","perm")} for f in "TE"}
    nrej={f:{k:[] for k in ("ind","cam","perm")} for f in "TE"}
    for s in range(NS):
        cond=np.concatenate([rng.permutation(6) for _ in range(4)])
        beta,XtXi,s2,Rz=fit(YY,cond)
        for k in range(1,6):
            z=zstat(beta,XtXi,s2,k); pi,pc,pp=pvals(z,Rz)
            for fam,mask in (("T",isT),("E",~isT)):
                for key,pv in (("ind",pi),("cam",pc),("perm",pp)):
                    v=pv[mask]; allp[fam][key].append(v); nrej[fam][key].append(int((false_discovery_control(v,method="bh")<.05).sum()))
    print(lab)
    for fam,l in (("T","task sets"),("E","essential sets")):
        for key,l2 in (("ind","independence"),("cam","CAMERA-like"),("perm","gene-label perm")):
            v=np.concatenate(allp[fam][key]); print("   ",l,l2,"FPR@.05=%.3f meanBH=%.2f"%(np.mean(v<.05),np.mean(nrej[fam][key])),flush=True)
