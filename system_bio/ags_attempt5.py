"""E14 hub-preserving null (relabel drugs): is the consensus better than chance when positives keep their drug-hub structure?  E15 skill beyond the PI hub for the best consensus."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
rk = rankdata
R3 = [ex(f, "rand") - ex(f, "ss") for f in FOLD]
V3 = [vf(f, "ss") for f in FOLD]
X3 = [ex(f, "ss") for f in FOLD]
C = {"vote fraction, link-only": V3[0], "vote fraction, topology-only": V3[1], "vote fraction, topo-and-link": V3[2],
     "paper rule, link-only": R3[0], "paper rule, topology-only": R3[1], "paper rule, topo-and-link": R3[2],
     "calibrated excess, wrong sign (higher=syn)": X3[0],
     "rank-mean vote, 3 strategies": np.mean([rk(v) for v in V3], 0),
     "rank-mean paper rule, 3 strategies": np.mean([rk(v) for v in R3], 0),
     "rank-mean paper rule 3 + vote 3": np.mean([rk(v) for v in (*R3, *V3)], 0)}
names = list(C)
RK = np.array([rk(C[n]) for n in names])
drugs = sorted({d for p in pairs for d in p})
pidx = {p: i for i, p in enumerate(pairs)}
obs_pairs = [tuple(p) for p in obs]
n_pos, n_neg = 6, 147
auc_from = lambda idx: (RK[:, idx].sum(1) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
obs_auc = auc_from(np.array([pidx[frozenset(p)] for p in obs_pairs]))
B = 50000
perm = np.empty((B, len(names)))
for b in range(B):
    m = dict(zip(drugs, rng.permutation(drugs)))
    perm[b] = auc_from(np.array([pidx[frozenset((m[a], m[c]))] for a, c in obs_pairs]))
mx = perm.max(1)
print("== E14  Hub-preserving null: relabel the 18 drugs at random (the 6 positives keep their shape: one drug in 4 of them)")
print(f"  {'score':44}{'AUC':>7}  {'p (hub null)':>13}  {'family-wise p':>14}")
for i in np.argsort(-obs_auc):
    print(f"  {names[i]:44}{obs_auc[i]:7.3f}  {np.mean(perm[:, i] >= obs_auc[i]):13.4f}  {np.mean(mx >= obs_auc[i]):14.4f}")
print(f"  (null AUC of a random score under this null: mean {perm.mean():.3f}; typical best-of-{len(names)} {np.percentile(mx, 95):.3f} at the 95th percentile)")

print("\n== E15  Skill beyond the PI hub for the best consensus (rank-mean paper rule 3 + vote 3)")
s = C["rank-mean paper rule 3 + vote 3"]
pi = np.array(["PI" in p for p in pairs])
for label, mask in (("within PI pairs (17 pairs, 4 positives)", pi), ("within other pairs (136 pairs, 2 positives)", ~pi)):
    yy, ss_ = y[mask], s[mask]
    b = roc_auc_score(yy, ss_)
    p = np.mean([roc_auc_score(rng.permutation(yy), ss_) >= b for _ in range(5000)])
    print(f"  {label:46} AUC {b:.2f}  permutation p={p:.3f}")
pos_rank = sorted(int(rk(-s, method='min')[pairs.index(p)]) for p in obs)
print(f"  ranks of the six synergies among 153: {pos_rank}")
