"""E23 leave-one-positive-out validation of the SELECTION procedure: pick the best of 16 label-free candidate scores on five positives, evaluate on the held-out sixth."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
exm = lambda f, k, m: load(f, k, 150, m)["Response excess over subset"].reindex(pairs).fillna(0).values
vfm = lambda f, k, m: frac(load(f, k, 150, m, "modelwise")).fillna(0).values
S = {}
for m in ("bliss", "hsa"):
    for f in FOLD:
        S[f"rule {m} {f}"] = exm(f, "rand", m) - exm(f, "ss", m)
        S[f"vote {m} {f}"] = vfm(f, "ss", m)
mr = lambda keys: np.mean([rk(S[k]) for k in keys], 0)
B6 = [k for k in S if "bliss" in k]
H6 = [k for k in S if "hsa" in k]
S["consensus Bliss"], S["consensus HSA"], S["consensus Bliss+HSA"], S["consensus Bliss rules + HSA votes"] = mr(B6), mr(H6), mr(B6 + H6), mr([k for k in B6 if k.startswith("rule")] + [k for k in H6 if k.startswith("vote")])
names = list(S)
pos, neg = np.where(y)[0], np.where(~y)[0]
pctl = lambda s, j: np.mean(s[neg] < s[j]) + .5 * np.mean(s[neg] == s[j])
picked, pct = [], []
for j in pos:
    rest = [k for k in pos if k != j]
    # AUC on the five remaining positives vs all negatives (negatives are not held out)
    best = max(names, key=lambda n: np.mean([pctl(S[n], k) for k in rest]))
    picked.append(best)
    pct.append(pctl(S[best], j))
print(f"== E23  Selection among {len(names)} label-free scores with leave-one-positive-out")
for j, b, p in zip(pos, picked, pct):
    print(f"  held-out {'~'.join(sorted(pairs[j])):6}: selected '{b}' -> held-out positive beats {p:.0%} of negatives")
print(f"  mean held-out percentile (= LOPO AUC of the selection procedure): {np.mean(pct):.3f}")
single = {n: np.mean([pctl(S[n], j) for j in pos]) for n in names}
print(f"  for comparison, median AUC of all {len(names)} candidates: {np.median(list(single.values())):.3f}; best {max(single.values()):.3f} ({max(single, key=single.get)}); worst {min(single.values()):.3f}")
