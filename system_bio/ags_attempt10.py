"""E21 where does the consensus signal come from (additive drug main effects vs pair-specific interaction)?  E22 learned weights vs equal weights (leave-one-positive-out)."""
import warnings
from scipy.stats import rankdata
from sklearn.linear_model import LogisticRegression
from ags_eval import *

warnings.filterwarnings("ignore")
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
exn = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
comp = {}
for f in FOLD:
    comp[f"rule {f}"] = rk(exn(f, "rand") - exn(f, "ss"))
    comp[f"vote {f}"] = rk(vfn(f, "ss"))
F = np.array(list(comp.values())).T  # 153 x 6 rank features
cons = F.mean(1)

drugs = sorted({d for p in pairs for d in p})
D = np.array([[d in p for d in drugs] for p in pairs], float)  # 153 x 18 incidence
print("== E21  Additive decomposition of each score: score(a,b) = drug_a + drug_b + pair-specific residual (fitted on the scores only, no labels)")
print(f"  {'score':26}{'AUC full':>10}{'AUC drug main effects':>24}{'AUC pair-specific residual':>28}{'share of variance from drugs':>30}")
for name, s in list(comp.items()) + [("consensus", cons)]:
    coef, *_ = np.linalg.lstsq(D, s, rcond=None)
    add = D @ coef
    res = s - add
    print(f"  {name:26}{auc(s):10.3f}{auc(add):24.3f}{auc(res):28.3f}{np.var(add) / np.var(s):30.2f}")
top_drug = sorted(zip(np.linalg.lstsq(D, cons, rcond=None)[0], drugs), reverse=True)[:4]
print("  largest drug main effects in the consensus: " + ", ".join(f"{d} ({c:+.0f})" for c, d in top_drug))

print("\n== E22  Learned weights vs equal weights: leave-one-positive-out, held-out positive's percentile among the 147 negatives")
pos, neg = np.where(y)[0], np.where(~y)[0]
for C in (0.01, 0.1, 1.0):
    pct = []
    for j in pos:
        idx = np.concatenate([[k for k in pos if k != j], neg])
        yy = np.concatenate([np.ones(len(pos) - 1), np.zeros(len(neg))])
        m = LogisticRegression(C=C, class_weight="balanced", max_iter=1000).fit((F[idx] - F.mean(0)) / F.std(0), yy)
        sc = ((F - F.mean(0)) / F.std(0)) @ m.coef_[0]
        pct.append(np.mean(sc[neg] < sc[j]) + .5 * np.mean(sc[neg] == sc[j]))
    print(f"  ridge logistic, C={C:<5}: LOPO AUC {np.mean(pct):.3f}  (per held-out positive: {', '.join(f'{p:.2f}' for p in pct)})")
eq = [np.mean(cons[neg] < cons[j]) + .5 * np.mean(cons[neg] == cons[j]) for j in pos]
print(f"  equal-weight consensus (nothing learned): AUC {np.mean(eq):.3f}  (per positive: {', '.join(f'{p:.2f}' for p in eq)})")
pc = np.linalg.svd(F - F.mean(0), full_matrices=False)
pc1 = (F - F.mean(0)) @ pc[2][0]
pc1 = pc1 if np.corrcoef(pc1, cons)[0, 1] > 0 else -pc1
print(f"  label-free first principal component of the six scores: AUC {auc(pc1):.3f} (loadings {np.round(np.abs(pc[2][0]), 2).tolist()})")
