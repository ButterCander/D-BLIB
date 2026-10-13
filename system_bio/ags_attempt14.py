"""E28 early-stopped vs max-generation models; E29 split-half by GA run; E30 run-level ensemble-size curve for the vote and for the consensus."""
import re
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
rng = np.random.default_rng(0)
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
auc = lambda s: roc_auc_score(y, s)
M = pd.read_csv(f"{A}/results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_model_predictions.tab", sep="\t", index_col=0)
P = [re.fullmatch(r"network_run_(\d+)__G(\d+)_M(\d+)", n) for n in M.index]
run, gen = np.array([int(p.group(1)) for p in P]), np.array([int(p.group(2)) for p in P])
M.columns = [pk(c) for c in M.columns]
M = M[pairs].values
vote = lambda m: np.nansum(m, 0) / np.maximum((~np.isnan(m)).sum(0), 1)
exn = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k: frac(load(f, k, 150, "bliss", "modelwise")).fillna(0).values
RR = [rk(exn(f, "rand") - exn(f, "ss")) for f in FOLD]
V_other = [rk(vfn(f, "ss")) for f in FOLD[1:]]
cons = lambda v: np.mean([*RR, *V_other, rk(v)], 0)

print("== E28  Models that reached the fitness target early (generation < 19) vs those that ran to the generation cap (G19)")
early, late = gen < 19, gen == 19
print(f"  early: {early.sum()} models, vote AUC {auc(vote(M[early])):.3f} | late: {late.sum()} models, vote AUC {auc(vote(M[late])):.3f}")
bal = [auc(vote(M[rng.choice(np.where(late)[0], early.sum(), replace=False)])) for _ in range(300)]
print(f"  late models, subsampled to {early.sum()}: AUC {np.mean(bal):.3f} +/- {np.std(bal):.3f}   (early models' AUC sits at the {np.mean(np.array(bal) < auc(vote(M[early]))):.0%} percentile of that distribution)")
print(f"  consensus with early-only votes {auc(cons(vote(M[early]))):.3f}, late-only votes {auc(cons(vote(M[late]))):.3f}, all {auc(cons(vote(M))):.3f}")

print("\n== E29  Split-half by GA run (75 runs vs 75 runs, 3 models each), 300 random splits")
ua = np.unique(run)
vh, ch = [], []
for _ in range(300):
    perm = rng.permutation(ua)
    for half in (perm[:75], perm[75:]):
        m = M[np.isin(run, half)]
        vh.append(auc(vote(m)))
        ch.append(auc(cons(vote(m))))
print(f"  vote AUC per half: mean {np.mean(vh):.3f}, sd {np.std(vh):.3f}, range {np.min(vh):.2f}-{np.max(vh):.2f}")
print(f"  consensus AUC per half: mean {np.mean(ch):.3f}, sd {np.std(ch):.3f}, range {np.min(ch):.2f}-{np.max(ch):.2f}")

print("\n== E30  Ensemble size in GA runs (each run contributes 3 models), draws without replacement")
print(f"  {'runs':>5}{'vote AUC':>18}{'consensus AUC':>18}")
for k in (5, 10, 25, 50, 75, 100, 125):
    v, c = [], []
    for _ in range(300):
        m = M[np.isin(run, rng.choice(ua, k, replace=False))]
        v.append(auc(vote(m)))
        c.append(auc(cons(vote(m))))
    print(f"  {k:5}{np.mean(v):10.3f} +/- {np.std(v):.3f}{np.mean(c):10.3f} +/- {np.std(c):.3f}")
print(f"  {150:5}{auc(vote(M)):10.3f}{'':8}{auc(cons(vote(M))):10.3f}")
