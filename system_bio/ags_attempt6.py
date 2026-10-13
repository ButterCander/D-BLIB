"""E16 ablation: is the consensus gain robust to which components are included (all 63 non-empty subsets of 3 paper-rule + 3 vote scores)?"""
import itertools
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
comp = {f"rule:{f}": rankdata(ex(f, "rand") - ex(f, "ss")) for f in FOLD}
comp.update({f"vote:{f}": rankdata(vf(f, "ss")) for f in FOLD})
keys = list(comp)
rows = []
for k in range(1, 7):
    for sub in itertools.combinations(keys, k):
        s = np.mean([comp[c] for c in sub], 0)
        rows.append((sub, roc_auc_score(y, s), average_precision_score(y, s), int((rankdata(-s, method="min")[y] <= 35).sum())))
a = np.array([r[1] for r in rows])
both = np.array([any(c.startswith("rule") for c in r[0]) and any(c.startswith("vote") for c in r[0]) for r in rows])
rule_only = np.array([all(c.startswith("rule") for c in r[0]) for r in rows])
vote_only = np.array([all(c.startswith("vote") for c in r[0]) for r in rows])
print(f"== E16  {len(rows)} subsets of the six component scores (rank-mean)")
for name, m in (("rule scores only", rule_only), ("vote scores only", vote_only), ("at least one rule AND one vote score", both)):
    print(f"  {name:40} n={m.sum():2}  AUC median {np.median(a[m]):.3f}  min {a[m].min():.3f}  max {a[m].max():.3f}   share with AUC >= 0.85: {np.mean(a[m] >= .85):.0%}")
print(f"  full six-score consensus: AUC {a[-1]:.3f}")
worst = sorted([r for r, b in zip(rows, both) if b], key=lambda r: r[1])[:3]
for r in worst:
    print(f"  weakest mixed subset: {', '.join(r[0])}  AUC {r[1]:.3f}")
