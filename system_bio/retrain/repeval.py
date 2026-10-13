"""Replication of the 32-node effect with NEW GA seeds, and random-node negative controls (link-only)."""
import os
import pickle
import warnings
import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score
from scipy.stats import rankdata

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
rk = rankdata
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]


def load(kind, runs):
    parts = [pickle.load(open(os.path.join(HERE, "out", f"balanced_{kind}_{r}.pkl"), "rb")) for r in runs]
    out = dict(parts[0])
    for k in ("fitness", "responses"):
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


pairs, e_r, v_r = scores(load("rand", (12, 24)))
y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
pos, neg = np.where(y)[0], np.where(~y)[0]


def row(label, kind, runs):
    if not all(os.path.exists(os.path.join(HERE, "out", f"balanced_{kind}_{r}.pkl")) for r in runs):
        return None
    d = load(kind, runs)
    _, e, v = scores(d)
    rule = e_r - e
    cons = rk(rule) + rk(v)
    print(f"{label:46} models {len(d['fitness']):3} | rule {roc_auc_score(y, rule):.3f} | vote {roc_auc_score(y, v):.3f} | consensus {roc_auc_score(y, cons):.3f} (PR {average_precision_score(y, cons):.3f}, top-35 {int((rk(-cons, method='min')[y] <= 35).sum())}/6)")
    return {"rule": rule, "vote": v, "cons": cons}


def row_tail(label, kind, skip_models):
    """the 36-run batch used GA seeds 36-71; seeds 36-47 duplicate runs 12-23 of the 24-run batch -> drop the first 12 runs (36 models)"""
    d = load(kind, (36,))
    d["responses"] = d["responses"][skip_models:]
    d["fitness"] = d["fitness"][skip_models:]
    _, e, v = scores(d)
    rule = e_r - e
    cons = rk(rule) + rk(v)
    print(f"{label:46} models {len(d['fitness']):3} | rule {roc_auc_score(y, rule):.3f} | vote {roc_auc_score(y, v):.3f} | consensus {roc_auc_score(y, cons):.3f} (PR {average_precision_score(y, cons):.3f}, top-35 {int((rk(-cons, method='min')[y] <= 35).sum())}/6)")
    return {"rule": rule, "vote": v, "cons": cons}


print("== A. replication with seeds 48-71 only (24 runs = 72 models each; the other 12 runs of the 36-run batch duplicated earlier seeds and are dropped)")
A = row_tail("24 nodes (baseline), seeds 48-71", "ss", 36)
B = row_tail("32 nodes (+8 expression nodes), seeds 48-71", "aug", 36)
if A and B:
    rng = np.random.default_rng(0)
    for k in ("rule", "vote", "cons"):
        d = []
        for _ in range(3000):
            i = np.concatenate([rng.choice(pos, 6), rng.choice(neg, 147)])
            d.append(roc_auc_score(y[i], B[k][i]) - roc_auc_score(y[i], A[k][i]))
        print(f"  difference (32 minus 24 nodes), {k}: {roc_auc_score(y, B[k]) - roc_auc_score(y, A[k]):+.3f}  95% [{np.percentile(d, 2.5):+.2f}, {np.percentile(d, 97.5):+.2f}]")
def pooled(label, kind):
    d1, d2, d3 = load(kind, (12,)), load(kind, (24,)), load(kind, (36,))
    d = dict(d1)
    d["responses"] = d1["responses"] + d2["responses"] + d3["responses"][36:]
    d["fitness"] = d1["fitness"] + d2["fitness"] + d3["fitness"][36:]
    _, e, v = scores(d)
    rule = e_r - e
    cons = rk(rule) + rk(v)
    print(f"{label:46} models {len(d['fitness']):3} | rule {roc_auc_score(y, rule):.3f} | vote {roc_auc_score(y, v):.3f} | consensus {roc_auc_score(y, cons):.3f} (PR {average_precision_score(y, cons):.3f}, top-35 {int((rk(-cons, method='min')[y] <= 35).sum())}/6)")


print("\n== B. pooled over all non-duplicated seeds (12 + 24 + 24 runs = 180 models each)")
pooled("24 nodes (baseline), all distinct seeds", "ss")
pooled("32 nodes, all distinct seeds", "aug")
print("\n== C. random-node negative controls vs real expression nodes (36 models each, GA seeds 12-23, paired with the first batch)")
row("24 nodes (baseline)", "ss", (12,))
row("32 nodes (+8 expression nodes)", "aug", (12,))
for c in ("ctl1", "ctl2", "ctl3"):
    row(f"24 + 8 RANDOM nodes with random states ({c})", c, (12,))
