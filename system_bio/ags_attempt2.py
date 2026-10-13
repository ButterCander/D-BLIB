"""E7 label-free calibration-specific component; E8 sensitivity of the paper's combination rule and leave-one-positive-out (LOPO) honesty check."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
pos, neg = np.where(y)[0], np.where(~y)[0]
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
z = lambda s: (s - s.mean()) / s.std()


def ci(s, B=2000, seed=2):
    r, o = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, len(pos)), r.choice(neg, len(neg))])
        o.append((roc_auc_score(y[i], s[i]), average_precision_score(y[i], s[i])))
    return np.percentile(o, [2.5, 97.5], axis=0)


print("== E7  Label-free 'calibration-specific' component: residual of calibrated score after regressing out the random-model score (higher = synergy)")
for f in FOLD:
    for name, g in (("ensemble excess", ex), ("vote fraction", vf)):
        s, r = g(f, "ss"), g(f, "rand")
        b = np.cov(s, r)[0, 1] / r.var(ddof=1)
        res = s - b * r
        lo, hi = ci(res)
        print(f"  {f:14} {name:16} slope {b:+.2f}  residual AUC {auc(res):.3f} [{lo[0]:.2f},{hi[0]:.2f}]   (calibrated alone {auc(s):.3f}, random alone {auc(r):.3f})")

print("\n== E8  The paper's combination score = calibrated + beta * random-model score (beta = -1 reported); evaluated as the paper does, lower = synergy -> score rand - ss")
TH = np.linspace(0, 2 * np.pi, 73)[:-1]
for f in FOLD:
    for name, g in (("ensemble excess", ex), ("vote fraction", vf)):
        s, r = g(f, "ss"), g(f, "rand")
        zs, zr = z(s), z(r)
        curve = np.array([auc(np.cos(t) * zs + np.sin(t) * zr) for t in TH])
        best = TH[curve.argmax()]
        # leave-one-positive-out: tune the angle on the other 5 positives (all negatives kept), score the held-out positive's percentile among negatives
        pct = []
        for j in pos:
            keep = np.array([k for k in range(len(y)) if k not in pos or k == j])
            rest = np.array([k for k in pos if k != j])
            idx = np.concatenate([rest, neg])
            yy = np.concatenate([np.ones(len(rest)), np.zeros(len(neg))])
            ac = [roc_auc_score(yy, (np.cos(t) * zs + np.sin(t) * zr)[idx]) for t in TH]
            t = TH[int(np.argmax(ac))]
            sc = np.cos(t) * zs + np.sin(t) * zr
            pct.append(np.mean(sc[neg] < sc[j]) + .5 * np.mean(sc[neg] == sc[j]))
        raw = r - s
        lo, hi = ci(raw)
        print(f"  {f:14} {name:16} rand-ss (paper rule) AUC {auc(raw):.3f} [{lo[0]:.2f},{hi[0]:.2f}] PR {average_precision_score(y, raw):.3f}"
              f" | best angle in-sample AUC {curve.max():.3f} | leave-one-positive-out AUC {np.mean(pct):.3f} | range over angles {curve.min():.2f}-{curve.max():.2f}")
    # top-35 hits under the paper rule (ensemble excess)
    raw = ex(f, "rand") - ex(f, "ss")
    print(f"  {f:14} top-35 pairs by paper rule contain {int((rankdata(-raw, method='min')[y] <= 35).sum())} of 6 synergies")
