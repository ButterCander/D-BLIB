"""Attempts to address the paper's limits using only the published result tables (no new simulation).
E1/E2 does calibration help? E3 missing-aware ensemble scores  E4 tie-breaking  E5 topology dependence  E6 skill beyond the PI hub."""
import warnings
from scipy.stats import rankdata, spearmanr
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
pos, neg = np.where(y)[0], np.where(~y)[0]
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
ap = lambda s: average_precision_score(y, s)
mk = lambda s, f=roc_auc_score: (lambda i: f(y[i], s[i]))


def dci(fa, fb, B=2000, seed=1):
    r, d = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, len(pos)), r.choice(neg, len(neg))])
        d.append(fa(i) - fb(i))
    return np.percentile(d, [2.5, 97.5])


def ci(s, B=2000, seed=2):
    r, o = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([r.choice(pos, len(pos)), r.choice(neg, len(neg))])
        o.append((roc_auc_score(y[i], s[i]), average_precision_score(y[i], s[i])))
    return np.percentile(o, [2.5, 97.5], axis=0)


print("== E1/E2  Does calibration help?  AUC (higher = synergy), Bliss, 150 runs; paired bootstrap CI of calibrated minus random-model score")
for f in FOLD:
    for name, g in (("ensemble excess", ex), ("vote fraction", vf)):
        s, r = g(f, "ss"), g(f, "rand")
        lo, hi = dci(mk(s), mk(r))
        print(f"  {f:14} {name:16} calibrated {auc(s):.3f}  random-model {auc(r):.3f}  diff {auc(s) - auc(r):+.3f} [{lo:+.2f},{hi:+.2f}]"
              f"   |  paper-style combo ss-rand: {auc(s - r):.3f} (higher=syn), {auc(r - s):.3f} (lower=syn)")

print("\n== E3  Missing-aware ensemble scores, link-only, 450 models (17% predictions missing)")
M = pd.read_csv(f"{A}/results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_model_predictions.tab", sep="\t", index_col=0)
M.columns = [pk(c) for c in M.columns]
M = M[pairs].values


def est(m, kind):
    ok, x = ~np.isnan(m), np.nan_to_num(m, nan=0.0)
    n, k = ok.sum(0), x.sum(0)
    if kind == "valid-only":
        return k / np.maximum(n, 1)
    if kind == "NaN as 0":
        return k / len(m)
    if kind == "shrunk (prior strength 20)":
        mu = k.sum() / n.sum()
        return (k + mu * 20) / (n + 20)
    if kind == "completeness-weighted":
        w = ok.mean(1)[:, None] * ok
        return (w * x).sum(0) / np.maximum(w.sum(0), 1e-9)
    if kind == "complete models only (>=80% valid)":
        return est(m[ok.mean(1) >= .8], "valid-only")


KINDS = ["valid-only", "NaN as 0", "shrunk (prior strength 20)", "completeness-weighted", "complete models only (>=80% valid)"]
print(f"  models with >=80% valid predictions: {(~np.isnan(M)).mean(1).__ge__(.8).sum()} of {len(M)};  pairs with fewer than 300 valid predictions: {((~np.isnan(M)).sum(0) < 300).sum()}")
print(f"  {'estimator':38}{'AUC [95% CI]':>22}{'PR AUC':>9}{'split-half reliability (Spearman-Brown)':>42}")
V = {}
for kd in KINDS:
    s = est(M, kd)
    V[kd] = s
    rho = []
    for _ in range(200):
        p = rng.permutation(len(M))
        r_ = spearmanr(est(M[p[:len(M) // 2]], kd), est(M[p[len(M) // 2:]], kd))[0]
        rho.append(2 * r_ / (1 + r_))
    lo, hi = ci(s)
    print(f"  {kd:38}{auc(s):8.3f} [{lo[0]:.2f},{hi[0]:.2f}]{ap(s):9.3f}{np.mean(rho):28.3f} +/- {np.std(rho):.3f}")

print("\n== E4  Tie-breaking the 105 zero-vote pairs and consensus (link-only); paired bootstrap vs plain vote fraction")
v, e = V["shrunk (prior strength 20)"], ex("link-only", "ss")
base = V["valid-only"]
for name, s in (("vote fraction (valid-only)", base), ("shrunk vote fraction", v), ("shrunk vote, ties broken by excess", v + 1e-6 * rankdata(e)),
                ("rank-mean(vote, excess)", rankdata(v) + rankdata(e))):
    lo, hi = ci(s)
    d1, d2 = dci(mk(s), mk(base)), dci(mk(s, average_precision_score), mk(base, average_precision_score))
    print(f"  {name:36} AUC {auc(s):.3f} [{lo[0]:.2f},{hi[0]:.2f}]  PR AUC {ap(s):.3f} [{lo[1]:.2f},{hi[1]:.2f}]   vs plain: dAUC [{d1[0]:+.2f},{d1[1]:+.2f}]  dPR [{d2[0]:+.2f},{d2[1]:+.2f}]")

print("\n== E5  How topology-driven are the predictions?  Spearman correlation across the 153 pairs")
VS = {f: vf(f, "ss") for f in FOLD}
ES = {f: ex(f, "ss") for f in FOLD}
for f in FOLD:
    print(f"  {f:14} calibrated vs random-model: vote fraction rho={spearmanr(vf(f, 'ss'), vf(f, 'rand'))[0]:.2f}, ensemble excess rho={spearmanr(ex(f, 'ss'), ex(f, 'rand'))[0]:.2f}")
for a_, b_ in (("link-only", "topology-only"), ("link-only", "topo-and-link"), ("topology-only", "topo-and-link")):
    print(f"  {a_} vs {b_}: vote fraction rho={spearmanr(VS[a_], VS[b_])[0]:.2f}, ensemble excess rho={spearmanr(ES[a_], ES[b_])[0]:.2f}")

print("\n== E6  Skill beyond the PI hub: AUC within PI-containing pairs (17 pairs, 4 positives) and within the other 136 pairs (2 positives); permutation p")
pi = np.array(["PI" in p for p in pairs])


def within(s, mask, B=4000):
    yy, ss_ = y[mask], s[mask]
    b = roc_auc_score(yy, ss_)
    return b, np.mean([roc_auc_score(rng.permutation(yy), ss_) >= b for _ in range(B)])


for name, s in (("ensemble excess, calibrated (link-only)", ex("link-only", "ss")), ("vote fraction, calibrated (link-only)", vf("link-only", "ss")),
                ("vote fraction, calibrated (topology-only)", vf("topology-only", "ss")), ("vote fraction, calibrated (topo-and-link)", vf("topo-and-link", "ss")),
                ("ensemble excess, random models (link-only)", ex("link-only", "rand")), ("vote fraction, random models (link-only)", vf("link-only", "rand"))):
    (a1, p1), (a2, p2) = within(s, pi), within(s, ~pi)
    print(f"  {name:44} within PI pairs AUC {a1:.2f} (p={p1:.2f})   within non-PI pairs AUC {a2:.2f} (p={p2:.3f})")
