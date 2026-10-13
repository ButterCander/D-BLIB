"""One-command reproduction of the headline numbers (reads only ags/ result tables). Run: python ags_headline.py"""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
rk = rankdata
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s, yy=y: roc_auc_score(yy, s)
top = lambda s, k=35: int((rk(-s, method="min")[y] <= k).sum())
pos, neg = np.where(y)[0], np.where(~y)[0]


def pair_ci(s, B=2000):
    r, o = np.random.default_rng(1), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, 6), r.choice(neg, 147)])
        o.append(roc_auc_score(y[i], s[i]))
    return np.percentile(o, [2.5, 97.5])


R = {f: ex(f, "rand") - ex(f, "ss") for f in FOLD}   # paper rule: random minus calibrated; NEGATIVE excess = synergy, so higher = synergy
V = {f: vf(f, "ss") for f in FOLD}
cons = np.mean([rk(v) for v in (*R.values(), *V.values())], 0)
print("1. paper rule (AUC, top-35 hits): " + "; ".join(f"{f} {auc(R[f]):.3f}, {top(R[f])}/6" for f in FOLD))
print("2. vote fraction AUC: " + "; ".join(f"{f} {auc(V[f]):.3f}" for f in FOLD))
print(f"3. consensus of 6 scores: AUC {auc(cons):.3f}, pair-bootstrap 95% {pair_ci(cons)[0]:.2f}-{pair_ci(cons)[1]:.2f}, PR AUC {average_precision_score(y, cons):.3f}, top-35 {top(cons)}/6, top-50 {top(cons, 50)}/6")
print(f"4. trivial rule 'pair contains PI': AUC {auc(np.array(['PI' in p for p in pairs], float)):.3f}")
pi = np.array(["PI" in p for p in pairs])
print(f"5. consensus within PI pairs AUC {auc(cons[pi], y[pi]):.2f}; within other pairs {auc(cons[~pi], y[~pi]):.2f}")
drugs = sorted({d for p in pairs for d in p})
jk = [auc(cons[np.array([d not in p for p in pairs])], y[np.array([d not in p for p in pairs])]) for d in drugs]
print(f"6. drug jackknife of the consensus: min {min(jk):.3f}, median {np.median(jk):.3f}")
pidx = {p: i for i, p in enumerate(pairs)}
obs_pairs = [tuple(p) for p in obs]
perm_max, perm_own = [], []
fam = [R[f] for f in FOLD] + [V[f] for f in FOLD] + [cons]
RK = np.array([rk(s) for s in fam])
for _ in range(20000):
    m = dict(zip(drugs, rng.permutation(drugs)))
    idx = np.array([pidx[frozenset((m[a], m[c]))] for a, c in obs_pairs])
    a = (RK[:, idx].sum(1) - 21) / (6 * 147)
    perm_max.append(a.max()), perm_own.append(a[-1])
print(f"7. hub-preserving permutation: own p {np.mean(np.array(perm_own) >= auc(cons)):.4f}, best-of-7 family-wise p {np.mean(np.array(perm_max) >= auc(cons)):.4f}")
