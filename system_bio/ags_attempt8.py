"""E17 does the consensus survive adding the independent HSA synergy metric (hub-preserving permutation, bigger family)?  E18 which pairs would be prioritised, and how stable are their ranks?"""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
FOLD = ["link-only", "topology-only", "topo-and-link"]
exm = lambda f, k, m: load(f, k, 150, m)["Response excess over subset"].reindex(pairs).fillna(0).values
vfm = lambda f, k, m: frac(load(f, k, 150, m, "modelwise")).fillna(0).values
rk = rankdata
auc = lambda s: roc_auc_score(y, s)
S = {}
for m in ("bliss", "hsa"):
    for f in FOLD:
        S[f"rule {m} {f}"] = exm(f, "rand", m) - exm(f, "ss", m)
        S[f"vote {m} {f}"] = vfm(f, "ss", m)
print("== E17  AUC of each component and of consensus scores built from Bliss only, HSA only, and both")
for k, v in S.items():
    print(f"  {k:36} AUC {auc(v):.3f}")
mean_rank = lambda keys: np.mean([rk(S[k]) for k in keys], 0)
B6 = [k for k in S if "bliss" in k]
H6 = [k for k in S if "hsa" in k]
C = {"consensus Bliss (3 rule + 3 vote)": mean_rank(B6), "consensus HSA (3 rule + 3 vote)": mean_rank(H6), "consensus Bliss + HSA (12 scores)": mean_rank(B6 + H6),
     "consensus Bliss rules + HSA votes": mean_rank([k for k in B6 if k.startswith("rule")] + [k for k in H6 if k.startswith("vote")])}
for k, v in C.items():
    print(f"  {k:36} AUC {auc(v):.3f}  PR AUC {average_precision_score(y, v):.3f}  synergies in top 35: {int((rk(-v, method='min')[y] <= 35).sum())} of 6")

# hub-preserving permutation family test over all single components + consensus scores
names = list(S) + list(C)
RK = np.array([rk(S[n]) if n in S else rk(C[n]) for n in names])
drugs = sorted({d for p in pairs for d in p})
pidx = {p: i for i, p in enumerate(pairs)}
obs_pairs = [tuple(p) for p in obs]
a_from = lambda idx: (RK[:, idx].sum(1) - 21) / (6 * 147)
obs_auc = a_from(np.array([pidx[frozenset(p)] for p in obs_pairs]))
perm = np.array([a_from(np.array([pidx[frozenset((m_[a], m_[c]))] for a, c in obs_pairs]))
                 for m_ in (dict(zip(drugs, rng.permutation(drugs))) for _ in range(20000))])
mx = perm.max(1)
print(f"\n  hub-preserving permutation, family of {len(names)} scores (every component and consensus tried):")
for n in (list(C)):
    i = names.index(n)
    print(f"  {n:36} AUC {obs_auc[i]:.3f}  p (own) {np.mean(perm[:, i] >= obs_auc[i]):.4f}  family-wise p {np.mean(mx >= obs_auc[i]):.4f}")

print("\n== E18  Top-ranked pairs by the Bliss consensus and rank stability under resampling the 450 link-only models (vote part) ")
M = pd.read_csv(f"{A}/results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_model_predictions.tab", sep="\t", index_col=0)
M.columns = [pk(c) for c in M.columns]
M = M[pairs].values
vote = lambda m: np.nansum(m, 0) / np.maximum((~np.isnan(m)).sum(0), 1)
fixed = [rk(S[k]) for k in B6 if k != "vote bliss link-only"]
cons = lambda v: np.mean(fixed + [rk(v)], 0)
base = cons(vote(M))
top = np.argsort(-base)[:15]
R = np.array([rk(-cons(vote(M[rng.choice(len(M), len(M))])), method="min") for _ in range(300)])
print(f"  {'pair':8}{'observed':>9}{'rank':>6}{'rank 95% range':>16}")
for i in top:
    lo, hi = np.percentile(R[:, i], [2.5, 97.5])
    print(f"  {'~'.join(sorted(pairs[i])):8}{'yes' if y[i] else 'no':>9}{int(rk(-base, method='min')[i]):6}{f'{lo:.0f}-{hi:.0f}':>16}")
print(f"  observed synergies among the top 15: {int(y[top].sum())} of 6; non-observed pairs in the top 15 (candidate new combinations): {int((~y[top]).sum())}")
