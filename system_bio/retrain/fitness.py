"""Does calibration fitness relate to synergy-prediction skill? (tests the earlier hypothesis that early-stopped / well-fitting models carry no signal)"""
import os
import pickle
import sys
import warnings
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = sys.argv[1] if len(sys.argv) > 1 else "12"
OBS = {frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()}


def votes(d):
    """models x pairs matrix: 1 synergy, 0 not, nan missing (negative Bliss excess = synergy)"""
    lo, hi = d["min_output"], d["max_output"]
    n = lambda v: (v - lo) / (hi - lo)
    pairs = sorted({k for r in d["responses"] for k in r if "-" in k})
    M = np.full((len(d["responses"]), len(pairs)), np.nan)
    for i, r in enumerate(d["responses"]):
        for j, p in enumerate(pairs):
            a, b = p.split("-")
            if any(r.get(k) is None for k in (a, b, p)):
                continue
            M[i, j] = float(n(r[p]) - n(r[a]) * n(r[b]) < 0)
    return pairs, M


for mt, label in (("balanced", "link-only"), ("topology", "topology-only"), ("mixed", "topo-and-link")):
    d = pickle.load(open(os.path.join(HERE, "out", f"{mt}_ss_{RUNS}.pkl"), "rb"))
    pairs, M = votes(d)
    y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
    fit = np.array(d["fitness"])
    vote = lambda m: np.nansum(m, 0) / np.maximum((~np.isnan(m)).sum(0), 1)
    single = np.array([roc_auc_score(y, np.nan_to_num(r)) for r in M])
    rho = spearmanr(fit, single)[0]
    order = np.argsort(-fit)
    half = len(fit) // 2
    hi_fit, lo_fit = order[:half], order[half:]
    print(f"{label:14} models {len(fit)}, fitness mean {fit.mean():.3f} (min {fit.min():.3f}); "
          f"Spearman(fitness, single-model AUC) = {rho:+.2f}; ensemble vote AUC all {roc_auc_score(y, vote(M)):.3f} | "
          f"fittest half {roc_auc_score(y, vote(M[hi_fit])):.3f} | least-fit half {roc_auc_score(y, vote(M[lo_fit])):.3f}; "
          f"fraction of fitness >= 0.99: {(fit >= .99).mean():.0%}")
