"""Calibration-data ablation (link-only, 12 GA runs = 36 models each): all 24 steady-state nodes vs two random subsets of 12 nodes."""
import os
import pickle
import warnings
import numpy as np
from sklearn.metrics import roc_auc_score
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]
ld = lambda kind: pickle.load(open(os.path.join(HERE, "out", f"balanced_{kind}_12.pkl"), "rb"))


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


pairs, e_r, v_r = scores(ld("rand"))
y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
rows = {}
for label, kind in (("all 24 nodes", "ss"), ("12 nodes, subset 1", "sub12_1"), ("12 nodes, subset 2", "sub12_2")):
    if not os.path.exists(os.path.join(HERE, "out", f"balanced_{kind}_12.pkl")):
        continue
    _, e, v = scores(ld(kind))
    rows[label] = (e_r - e, v)
    print(f"{label:20} paper rule AUC {roc_auc_score(y, e_r - e):.3f} | vote fraction AUC {roc_auc_score(y, v):.3f}")
names = list(rows)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        print(f"  rank correlation of vote fractions, {names[i]} vs {names[j]}: {spearmanr(rows[names[i]][1], rows[names[j]][1])[0]:.2f}; paper-rule scores: {spearmanr(rows[names[i]][0], rows[names[j]][0])[0]:.2f}")
