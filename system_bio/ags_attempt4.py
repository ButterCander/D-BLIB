"""E12 decide the score-sign convention with the paper's CASCADE 1.0 result (AUC 0.97); E13 family-wise permutation test for the best of all label-free scores tried."""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)

print("== E12  CASCADE 1.0 (21 pairs, 4 observed synergies; paper reports ROC AUC 0.97, PR AUC 0.91), link-only 50 runs")


def load1(kind, metric, level):
    d = pd.read_csv(f"{A}/results/link-only/cascade_1.0_{kind}_50sim_fixpoints_{metric}_{level}_synergies.tab", sep="\t", index_col=0)
    d.index = [pk(i) for i in d.index]
    return d


obs1 = {frozenset(l.strip().split("~")) for l in open(f"{A}/data/observed_synergies_cascade_1.0") if l.strip()}
for metric in ("bliss", "hsa"):
    mw = load1("ss", metric, "modelwise")
    p1 = list(mw.index)
    y1 = np.array([p in obs1 for p in p1])
    e_ss = load1("ss", metric, "ensemblewise").iloc[:, 0].reindex(p1).fillna(0).values
    e_r = load1("rand", metric, "ensemblewise").iloc[:, 0].reindex(p1).fillna(0).values
    v_ss = (mw.Synergies / (mw.Synergies + mw["Non-synergies"])).fillna(0).values
    f = lambda s: f"AUC {roc_auc_score(y1, s):.2f} PR {average_precision_score(y1, s):.2f}"
    print(f"  {metric}: {len(p1)} pairs, {y1.sum()} positives | calibrated excess higher=syn: {f(e_ss)} | lower=syn: {f(-e_ss)} | rule rand-ss higher=syn: {f(e_r - e_ss)} | ss-rand higher=syn: {f(e_ss - e_r)} | vote fraction: {f(v_ss)}"
          f" | rank-mean(rule, vote): {f(rankdata(e_r - e_ss) + rankdata(v_ss))}")

print("\n== E13  Family-wise permutation test: is the best of all label-free scores tried better than chance after accounting for the search?")
FOLD = ["link-only", "topology-only", "topo-and-link"]
ex = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vf = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
rk = rankdata
C = {}
for f in FOLD:
    C[f"vote fraction, {f}"] = vf(f, "ss")
    C[f"calibrated excess raw, {f}"] = ex(f, "ss")
    C[f"paper rule rand-ss, {f}"] = ex(f, "rand") - ex(f, "ss")
R3 = [ex(f, "rand") - ex(f, "ss") for f in FOLD]
V3 = [vf(f, "ss") for f in FOLD]
X3 = [ex(f, "ss") for f in FOLD]
C["mean vote, 3 strategies"] = np.mean(V3, 0)
C["rank-mean vote, 3 strategies"] = np.mean([rk(v) for v in V3], 0)
C["rank-mean paper rule, 3 strategies"] = np.mean([rk(v) for v in R3], 0)
C["rank-mean paper rule 3 + vote 3"] = np.mean([rk(v) for v in (*R3, *V3)], 0)
C["rank-mean vote 3 + raw excess 3"] = np.mean([rk(v) for v in (*V3, *X3)], 0)
C["rank-mean(vote, excess), link-only"] = rk(V3[0]) + rk(X3[0])
names = list(C)
RK = np.array([rk(C[n]) for n in names])
n_pos, n_neg = int(y.sum()), int((~y).sum())
obs_auc = np.array([(RK[i][y].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg) for i in range(len(names))])
B = 50000
perm_auc = np.empty((B, len(names)))
for b in range(B):
    idx = rng.choice(len(y), n_pos, replace=False)
    perm_auc[b] = (RK[:, idx].sum(1) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
mx = perm_auc.max(1)
print(f"  {len(names)} candidate scores tried (all label-free).")
for i in np.argsort(-obs_auc):
    print(f"  {names[i]:40} AUC {obs_auc[i]:.3f}   single-test p={np.mean(perm_auc[:, i] >= obs_auc[i]):.4f}   family-wise p={np.mean(mx >= obs_auc[i]):.4f}")
