"""Sanity checks on the 6-positive benchmark: trivial 'contains PI' baseline, ranks of positives, stability over ensemble size, delta-AUC CI."""
from ags_eval import *

rng = np.random.default_rng(1)
auc = lambda s, yy=y: roc_auc_score(yy, s)
print("positives:", sorted("~".join(sorted(p)) for p in obs), "  (PI in", sum("PI" in p for p in obs), "of 6)")

pi = np.array(["PI" in p for p in pairs], float)
print(f"trivial baseline 'pair contains PI': AUC {auc(pi):.3f}, PR AUC {average_precision_score(y, pi):.3f}")

# drug-level additive baseline built only from model outputs (no labels): mean vote fraction of each drug's pairs
vf = frac(load("link-only", "ss", 150, level="modelwise")).fillna(0)
drugs = sorted({d for p in pairs for d in p})
m = {d: np.mean([v for p, v in vf.items() if d in p]) for d in drugs}
add = np.array([sum(m[d] for d in p) for p in pairs])
print(f"additive single-drug baseline from model outputs: AUC {auc(add):.3f}")

print("\nvote fraction (ss) per setup: positives' rank among 153 (1=top) and #pairs with score 0")
for folder in ("link-only", "topology-only", "topo-and-link"):
    s = frac(load(folder, "ss", 150, level="modelwise")).fillna(0).values
    rk = (-s).argsort().argsort() + 1
    print(f"  {folder:14} ranks of positives {sorted(rk[y])}   zeros: {(s == 0).sum()}   AUC {auc(s):.3f}")

print("\nstability over ensemble size (link-only, AUC; ensemble excess | vote fraction)")
for n in (50, 100, 150, 200):
    e = load("link-only", "ss", n)["Response excess over subset"].reindex(pairs).fillna(-9).values
    v = frac(load("link-only", "ss", n, level="modelwise")).fillna(0).values
    print(f"  {n:3} runs ({3 * n} models): {auc(e):.3f} | {auc(v):.3f}")

e = load("link-only", "ss", 150)["Response excess over subset"].reindex(pairs).fillna(-9).values
pos, neg = np.where(y)[0], np.where(~y)[0]
d = []
for _ in range(5000):
    i = np.concatenate([rng.choice(pos, 6), rng.choice(neg, len(neg))])
    d.append(auc(v[i], y[i]) - auc(e[i], y[i]))
print(f"\nvote fraction minus ensemble excess (link-only): dAUC {auc(v) - auc(e):+.3f}, 95% CI [{np.percentile(d, 2.5):+.2f},{np.percentile(d, 97.5):+.2f}]")
