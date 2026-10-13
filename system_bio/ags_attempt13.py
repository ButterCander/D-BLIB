"""E25 aggregation/weight sensitivity of the consensus; E26 drug-cluster bootstrap CI; E27 negative control (random-model votes instead of calibrated votes)."""
import warnings
from scipy.stats import rankdata, gmean
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
exn = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
hits = lambda s, k=35: int((rk(-s, method="min")[y] <= k).sum())
R = [exn(f, "rand") - exn(f, "ss") for f in FOLD]
V = [vfn(f, "ss") for f in FOLD]
VR = [vfn(f, "rand") for f in FOLD]
RR, RV = [rk(v) for v in R], [rk(v) for v in V]
z = lambda v: (v - v.mean()) / v.std()

print("== E25a  Aggregation method (same six scores)")
allr = np.array([*RR, *RV])
for name, s in (("mean of ranks (used)", allr.mean(0)), ("median of ranks", np.median(allr, 0)), ("geometric mean of ranks", gmean(allr, axis=0)),
                ("mean of z-scores", np.mean([z(v) for v in (*R, *V)], 0)), ("best (max) rank only", allr.max(0)), ("worst (min) rank only", allr.min(0))):
    print(f"  {name:28} AUC {auc(s):.3f}  PR {average_precision_score(y, s):.3f}  top-35 hits {hits(s)}")

print("\n== E25b  Weight on the rule block vs the vote block: score = (1-w) * mean rank(vote x3) + w * mean rank(rule x3)")
print("  w:   " + "".join(f"{w:6.1f}" for w in np.linspace(0, 1, 11)))
print("  AUC: " + "".join(f"{auc((1 - w) * np.mean(RV, 0) + w * np.mean(RR, 0)):6.3f}" for w in np.linspace(0, 1, 11)))

print("\n== E26  Drug-cluster bootstrap (resample the 18 drugs, keep the pairs among resampled drugs): 95% interval of AUC")
drugs = sorted({d for p in pairs for d in p})
pidx = {p: i for i, p in enumerate(pairs)}
cons = allr.mean(0)
S = {"paper rule, parameters only": R[0], "paper rule, topology + parameters": R[2], "vote fraction, parameters only": V[0], "consensus of 6 scores": cons}
out = {k: [] for k in S}
for _ in range(2000):
    d = rng.choice(drugs, len(drugs))
    idx = [pidx[frozenset((d[i], d[j]))] for i in range(len(d)) for j in range(i + 1, len(d)) if d[i] != d[j]]
    idx = np.array(idx)
    if y[idx].sum() == 0 or y[idx].sum() == len(idx):
        continue
    for k, s in S.items():
        out[k].append(roc_auc_score(y[idx], s[idx]))
for k, v in out.items():
    print(f"  {k:36} point {auc(S[k]):.3f}  cluster-bootstrap median {np.median(v):.3f}  95% [{np.percentile(v, 2.5):.2f}, {np.percentile(v, 97.5):.2f}]  ({len(v)} valid resamples)")

print("\n== E27  Negative control: replace the calibrated vote fractions by random-model vote fractions")
for name, s in (("consensus (3 rule + 3 calibrated votes)", cons), ("3 rule + 3 RANDOM-model votes", np.mean([*RR, *[rk(v) for v in VR]], 0)),
                ("3 calibrated votes only", np.mean(RV, 0)), ("3 random-model votes only", np.mean([rk(v) for v in VR], 0)), ("3 rule only", np.mean(RR, 0))):
    print(f"  {name:42} AUC {auc(s):.3f}  PR {average_precision_score(y, s):.3f}  top-35 hits {hits(s)}")
