"""Does adding expression-derived calibration nodes (24 -> 32) change link-only predictions? Same GA seeds, paired comparison."""
import os
import pickle
import sys
import warnings
import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else "12,24"
rk = rankdata
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]


def load(kind):
    parts = [pickle.load(open(os.path.join(HERE, "out", f"balanced_{kind}_{r}.pkl"), "rb")) for r in RUNS.split(",")]
    out = dict(parts[0])
    for k in ("fitness", "names", "responses"):
        out[k] = [x for p in parts for x in p[k]]
    return out


def scores(d):
    lo, hi = d["min_output"], d["max_output"]
    n = lambda v: (v - lo) / (hi - lo)
    resp = d["responses"]
    pairs = sorted({k for r in resp for k in r if "-" in k})
    ens, vote = [], []
    for p in pairs:
        a, b = p.split("-")
        m = lambda k: np.mean([r[k] for r in resp if r.get(k) is not None] or [np.nan])
        ens.append(n(m(p)) - n(m(a)) * n(m(b)))
        s = [n(r[p]) - n(r[a]) * n(r[b]) < 0 for r in resp if all(r.get(k) is not None for k in (a, b, p))]
        vote.append(np.mean(s) if s else 0.0)
    return pairs, np.nan_to_num(np.array(ens)), np.array(vote)


pairs, e_r, v_r = scores(load("rand"))
y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
pos, neg = np.where(y)[0], np.where(~y)[0]
S = {}
for label, kind in (("24 nodes (baseline)", "ss"), ("32 nodes (+8 from AGS expression)", "aug")):
    d = load(kind)
    _, e, v = scores(d)
    rule = e_r - e
    S[label] = {"rule": rule, "vote": v, "cons": rk(rule) + rk(v)}
    print(f"{label:34} models {len(d['fitness'])}, fitness {np.mean(d['fitness']):.3f} | rule AUC {roc_auc_score(y, rule):.3f} | vote AUC {roc_auc_score(y, v):.3f} | "
          f"2-score consensus AUC {roc_auc_score(y, S[label]['cons']):.3f} PR {average_precision_score(y, S[label]['cons']):.3f}, top-35 hits {int((rk(-S[label]['cons'], method='min')[y] <= 35).sum())}/6")
a, b = list(S.values())
print(f"rank correlation between the two ensembles: rule {spearmanr(a['rule'], b['rule'])[0]:.2f}, vote {spearmanr(a['vote'], b['vote'])[0]:.2f}")
rng = np.random.default_rng(0)
for k in ("rule", "vote", "cons"):
    d = []
    for _ in range(3000):
        i = np.concatenate([rng.choice(pos, 6), rng.choice(neg, 147)])
        d.append(roc_auc_score(y[i], b[k][i]) - roc_auc_score(y[i], a[k][i]))
    print(f"  AUC difference (32 minus 24 nodes), {k}: {roc_auc_score(y, b[k]) - roc_auc_score(y, a[k]):+.3f}  95% [{np.percentile(d, 2.5):+.2f}, {np.percentile(d, 97.5):+.2f}]")
