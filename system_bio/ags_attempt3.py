"""E9 stability of the paper rule over ensemble size; E10 model-wise weighted rule (1-w)*p_cal + w*p_rand; E11 label-free consensus under the paper rule."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
pos, neg = np.where(y)[0], np.where(~y)[0]
FOLD = ["link-only", "topology-only", "topo-and-link"]
auc = lambda s: roc_auc_score(y, s)
ap = lambda s: average_precision_score(y, s)
exn = lambda f, k, n: load(f, k, n, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k, n: frac(load(f, k, n, "bliss", "modelwise")).fillna(0).values


def ci(s, B=2000, seed=2):
    r, o = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, len(pos)), r.choice(neg, len(neg))])
        o.append((roc_auc_score(y[i], s[i]), average_precision_score(y[i], s[i])))
    return np.percentile(o, [2.5, 97.5], axis=0)


print("== E9  Paper rule (random-model minus calibrated ensemble excess; higher = synergy) vs number of evolution runs")
for f in FOLD:
    row = []
    for n in (50, 150):
        try:
            s = exn(f, "rand", n) - exn(f, "ss", n)
            row.append(f"{n} runs: AUC {auc(s):.3f} PR {ap(s):.3f}")
        except FileNotFoundError:
            row.append(f"{n} runs: no random-model file")
    print(f"  {f:14} " + "   ".join(row))
for n in (50, 100, 150, 200):
    s = exn("link-only", "ss", n)
    print(f"  link-only calibrated-only excess, {n} runs: AUC {auc(s):.3f}  (random-model file available for 50 and 150 only)")

print("\n== E10  Model-wise weighted rule (1-w)*p_calibrated + w*p_random, AUC (higher = synergy)")
ws = np.linspace(0, 1, 11)
for f in FOLD:
    pc, pr = vfn(f, "ss", 150), vfn(f, "rand", 150)
    cur = [auc((1 - w) * pc + w * pr) for w in ws]
    print(f"  {f:14} w=0: {cur[0]:.3f}  w=0.5: {cur[5]:.3f}  w=1: {cur[10]:.3f}   min {min(cur):.3f} max {max(cur):.3f}")

print("\n== E11  Label-free consensus of the three calibration strategies under the paper rule (rank-mean of random-minus-calibrated excess)")
R = {f: exn(f, "rand", 150) - exn(f, "ss", 150) for f in FOLD}
for name, s in (("link-only alone", R["link-only"]), ("topology-only alone", R["topology-only"]), ("topo-and-link alone", R["topo-and-link"]),
                ("rank-mean of the three strategies", np.mean([rankdata(v) for v in R.values()], 0)),
                ("rank-mean of 3 strategies + 3 vote fractions", np.mean([rankdata(v) for v in (*R.values(), *[vfn(f, 'ss', 150) for f in FOLD])], 0))):
    lo, hi = ci(s)
    print(f"  {name:46} AUC {auc(s):.3f} [{lo[0]:.2f},{hi[0]:.2f}]  PR AUC {ap(s):.3f} [{lo[1]:.2f},{hi[1]:.2f}]  synergies in top 35: {int((rankdata(-s, method='min')[y] <= 35).sum())} of 6")
