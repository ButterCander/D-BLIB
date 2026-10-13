"""Evaluate the freshly re-simulated ensembles with the PRE-SPECIFIED recipes (nothing here was tuned on the new data):
paper rule = random-minus-calibrated ensemble Bliss excess; vote fraction of calibrated models; consensus = rank-mean of 3 rules + 3 votes."""
import os
import pickle
import sys
import warnings
import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else "24"
MT = {"link-only": "balanced", "topology-only": "topology", "topo-and-link": "mixed"}
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]
rk = rankdata


def load(mt, kind):
    """pool independent-seed batches, e.g. RUNS='12,24' -> 36 runs = 108 models per condition"""
    parts = [pickle.load(open(os.path.join(HERE, "out", f"{mt}_{kind}_{r}.pkl"), "rb")) for r in RUNS.split(",")]
    out = dict(parts[0])
    for k in ("fitness", "names", "responses"):
        out[k] = [x for p in parts for x in p[k]]
    return out


def bliss(r1, r2, rc, lo, hi):
    n = lambda v: (v - lo) / (hi - lo)
    return n(rc) - n(r1) * n(r2)


def scores(d):
    """ensemble-wise Bliss excess (mean responses) and model-wise vote fraction for every pair"""
    resp = d["responses"]
    lo, hi = d["min_output"], d["max_output"]
    pairs = sorted({k for r in resp for k in r if "-" in k})
    ens, vote, missing = {}, {}, {}
    for p in pairs:
        a, b = p.split("-")
        mean = lambda k: np.mean([r[k] for r in resp if isinstance(r.get(k), (int, float)) and r.get(k) is not None] or [np.nan])
        ens[p] = bliss(mean(a), mean(b), mean(p), lo, hi)
        syn = non = 0
        for r in resp:
            if any(r.get(k) is None for k in (a, b, p)):
                continue
            s = bliss(r[a], r[b], r[p], lo, hi)
            syn += s < 0
            non += s >= 0
        vote[p] = syn / (syn + non) if syn + non else 0.0
        missing[p] = sum(any(r.get(k) is None for k in (a, b, p)) for r in resp) / len(resp)
    return pairs, ens, vote, missing


pairs, y = None, None
R, V, info = {}, {}, {}
for name, mt in MT.items():
    ss, rd = load(mt, "ss"), load(mt, "rand")
    p1, e_ss, v_ss, m_ss = scores(ss)
    p2, e_rd, v_rd, m_rd = scores(rd)
    if pairs is None:
        pairs = p1
        y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
    # negative excess = synergy; paper rule = random - calibrated (higher = synergy)
    R[name] = np.nan_to_num(np.array([e_rd[p] - e_ss[p] for p in pairs]))
    V[name] = np.array([v_ss[p] for p in pairs])
    info[name] = (len(ss["fitness"]), np.mean(ss["fitness"]), np.mean([m_ss[p] for p in pairs]), np.mean(rd["fitness"]), V[name], np.array([v_rd[p] for p in pairs]))
auc = lambda s: roc_auc_score(y, s)
print(f"{len(pairs)} pairs, {y.sum()} observed synergies, {RUNS} GA runs per condition\n")
for name, (n, fit, miss, fit_r, v, vr) in info.items():
    print(f"{name:14} models {n}, mean fitness calibrated {fit:.3f} / random {fit_r:.3f}, mean missing prediction rate {miss:.1%}")
print("\npaper rule AUC:      " + "  ".join(f"{k} {auc(v):.3f}" for k, v in R.items()))
print("vote fraction AUC:   " + "  ".join(f"{k} {auc(v):.3f}" for k, v in V.items()))
print("random-model votes:  " + "  ".join(f"{k} {auc(info[k][5]):.3f}" for k in V))
cons = np.mean([rk(v) for v in (*R.values(), *V.values())], 0)
hit = lambda s, k: int((rk(-s, method="min")[y] <= k).sum())
print(f"\nCONSENSUS (3 rule + 3 vote, pre-specified): AUC {auc(cons):.3f}, PR AUC {average_precision_score(y, cons):.3f}, hits top-35 {hit(cons, 35)}/{y.sum()}, top-50 {hit(cons, 50)}/{y.sum()}")
print("ranks of the observed synergies:", sorted((int(rk(-cons, method='min')[i]), pairs[i]) for i in np.where(y)[0]))
pi = np.array(["PI" in p.split("-") for p in pairs])
print(f"within PI pairs AUC {roc_auc_score(y[pi], cons[pi]):.2f}; within other pairs {roc_auc_score(y[~pi], cons[~pi]):.2f}; trivial 'contains PI' {auc(pi.astype(float)):.3f}")
pub = {  # published-table results from the earlier analysis for comparison
    "rule": {"link-only": .686, "topology-only": .734, "topo-and-link": .823}, "vote": {"link-only": .831, "topology-only": .757, "topo-and-link": .804}, "consensus": .917}
print(f"published tables for comparison: rule {pub['rule']}, vote {pub['vote']}, consensus {pub['consensus']}")
print("rank correlation of my new consensus with the published-table consensus is computed in compare.py")
np.save(os.path.join(HERE, "out", f"consensus_{RUNS}.npy"), cons)
pickle.dump({"pairs": pairs, "y": y, "R": R, "V": V}, open(os.path.join(HERE, "out", f"scores_{RUNS}.pkl"), "wb"))
