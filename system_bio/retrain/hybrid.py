"""New method test: ONE ensemble trained on steady-state + proliferation (hybrid) vs the paper's two-ensemble contrast (random minus calibrated).
usage: python hybrid.py <runs of ss/rand to pool, e.g. 12 or 12,24>   (hybrid ensembles come from *_hyb_12.pkl)"""
import os
import pickle
import sys
import warnings
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else "12"
MT = {"link-only": "balanced", "topology-only": "topology", "topo-and-link": "mixed"}
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]
rk = rankdata


def load(mt, kind, runs):
    parts = [pickle.load(open(os.path.join(HERE, "out", f"{mt}_{kind}_{r}.pkl"), "rb")) for r in runs.split(",")]
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


y = None
S = {}
for name, mt in MT.items():
    pairs, e_ss, v_ss = scores(load(mt, "ss", RUNS))
    _, e_rd, v_rd = scores(load(mt, "rand", RUNS))
    if y is None:
        y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
    S[name] = {"rule": e_rd - e_ss, "vote": v_ss}
    hp = os.path.join(HERE, "out", f"{mt}_hyb_12.pkl")
    if os.path.exists(hp):
        _, e_h, v_h = scores(load(mt, "hyb", "12"))
        S[name]["hyb_excess"] = -e_h          # negative excess = synergy
        S[name]["hyb_vote"] = v_h
        S[name]["hyb_contrast"] = e_rd - e_h  # random minus hybrid, same form as the paper's rule
auc = lambda s: roc_auc_score(y, s)
print(f"ss/rand ensembles pooled from runs {RUNS}; hybrid ensembles 12 runs (36 models). 153 pairs, {y.sum()} positives\n")
print(f"{'strategy':14}" + "".join(f"{k:>14}" for k in ("rule", "vote", "hyb_excess", "hyb_vote", "hyb_contrast")))
for name in MT:
    print(f"{name:14}" + "".join(f"{auc(S[name][k]):14.3f}" if k in S[name] else f"{'n/a':>14}" for k in ("rule", "vote", "hyb_excess", "hyb_vote", "hyb_contrast")))
if all("hyb_excess" in S[n] for n in MT):
    base = np.mean([rk(S[n][k]) for n in MT for k in ("rule", "vote")], 0)
    hyb = np.mean([rk(S[n][k]) for n in MT for k in ("hyb_excess", "hyb_vote")], 0)
    both = np.mean([rk(S[n][k]) for n in MT for k in ("rule", "vote", "hyb_excess", "hyb_vote")], 0)
    for label, s in (("consensus of paper-rule + vote (6)", base), ("consensus of hybrid excess + vote (6)", hyb), ("consensus of all four kinds (12)", both)):
        print(f"{label:40} AUC {auc(s):.3f}  PR AUC {average_precision_score(y, s):.3f}  top-35 hits {int((rk(-s, method='min')[y] <= 35).sum())}/{y.sum()}")
