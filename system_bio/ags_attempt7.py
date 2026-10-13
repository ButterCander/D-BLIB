"""E2/E7 redone with the correct sign convention (ensemble excess: NEGATIVE = synergy, confirmed on CASCADE 1.0; vote fraction: HIGHER = synergy)."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
pos, neg = np.where(y)[0], np.where(~y)[0]
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
mk = lambda s: (lambda i: roc_auc_score(y[i], s[i]))


def dci(fa, fb, B=2000, seed=1):
    r, d = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, len(pos)), r.choice(neg, len(neg))])
        d.append(fa(i) - fb(i))
    return np.percentile(d, [2.5, 97.5])


print("== Ensemble excess, correct sign (lower = synergy):  AUC of score = -excess")
for f in FOLD:
    ss, rd = ex(f, "ss"), ex(f, "rand")
    rule = rd - ss
    b = np.cov(ss, rd)[0, 1] / rd.var(ddof=1)
    resid = -(ss - b * rd)
    l1, h1 = dci(mk(rule), mk(-ss))
    l2, h2 = dci(mk(rule), mk(-rd))
    print(f"  {f:14} calibrated alone {auc(-ss):.3f} | random-model alone {auc(-rd):.3f} | paper rule (random - calibrated) {auc(rule):.3f}"
          f" | rule minus calibrated-alone {auc(rule) - auc(-ss):+.2f} [{l1:+.2f},{h1:+.2f}] | rule minus random-alone {auc(rule) - auc(-rd):+.2f} [{l2:+.2f},{h2:+.2f}]"
          f" | label-free residual (slope {b:.2f}) {auc(resid):.3f}")
print("\n== Vote fraction (higher = synergy)")
for f in FOLD:
    pc, pr = vf(f, "ss"), vf(f, "rand")
    b = np.cov(pc, pr)[0, 1] / pr.var(ddof=1)
    l, h = dci(mk(pc), mk(pr))
    print(f"  {f:14} calibrated {auc(pc):.3f} | random-model {auc(pr):.3f} | diff {auc(pc) - auc(pr):+.3f} [{l:+.2f},{h:+.2f}] | label-free residual {auc(pc - b * pr):.3f}")

print("\n== Ranks (1 = top) of the six synergies under each score; scores = paper rule (random - calibrated excess), calibrated vote fraction, six-score consensus")
R3 = [ex(f, "rand") - ex(f, "ss") for f in FOLD]
V3 = [vf(f, "ss") for f in FOLD]
cons = np.mean([rankdata(v) for v in (*R3, *V3)], 0)
rank = lambda s, p: int(rankdata(-s, method="min")[pairs.index(p)])
for p in sorted(obs, key=lambda p: rank(cons, p)):
    print(f"  {'~'.join(sorted(p)):7} PI={'yes' if 'PI' in p else 'no ':3} paper rule (parameters only) {rank(R3[0], p):3}  paper rule (topology+parameters) {rank(R3[2], p):3}  vote fraction (parameters only) {rank(V3[0], p):3}  consensus {rank(cons, p):3}")
