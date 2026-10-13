"""E24 drug jackknife: drop every pair containing drug d, recompute the AUC of each score on the remaining pairs."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
exn = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
R = [exn(f, "rand") - exn(f, "ss") for f in FOLD]
V = [vfn(f, "ss") for f in FOLD]
S = {"paper rule, parameters only": R[0], "paper rule, topology + parameters": R[2], "vote fraction, parameters only": V[0],
     "consensus of 6 scores": np.mean([rk(v) for v in (*R, *V)], 0)}
drugs = sorted({d for p in pairs for d in p})
print("== E24  AUC after removing all pairs that contain drug d (positives left in brackets)")
print(f"  {'drug':5}{'pos left':>9}" + "".join(f"{k[:30]:>32}" for k in S))
rows = []
for d in drugs:
    keep = np.array([d not in p for p in pairs])
    if y[keep].sum() == 0 or y[keep].sum() == keep.sum():
        continue
    a = [roc_auc_score(y[keep], s[keep]) for s in S.values()]
    rows.append(a)
    print(f"  {d:5}{int(y[keep].sum()):9}" + "".join(f"{v:32.3f}" for v in a))
rows = np.array(rows)
print(f"  {'min':5}{'':9}" + "".join(f"{v:32.3f}" for v in rows.min(0)))
print(f"  {'median':5}{'':9}" + "".join(f"{v:32.3f}" for v in np.median(rows, 0)))
