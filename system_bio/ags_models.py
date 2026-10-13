"""Ensemble variability: resample the 450 individual models (link-only, ss) instead of the 153 pairs."""
from ags_eval import *
from scipy.stats import spearmanr

M = pd.read_csv(f"{A}/results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_model_predictions.tab", sep="\t", index_col=0)
M.columns = [pk(c) for c in M.columns]
M = M[pairs].values  # models x pairs, 1=synergy 0=no, nan=no prediction
print(f"{M.shape[0]} models, {np.isnan(M).mean():.0%} predictions missing, models predicting >=1 synergy: {(np.nansum(M, 1) > 0).mean():.0%}")

vote = lambda m: np.nansum(m, 0) / np.maximum((~np.isnan(m)).sum(0), 1)
rng = np.random.default_rng(0)
full = vote(M)
print(f"full ensemble vote fraction: AUC {roc_auc_score(y, full):.3f}, pairs with 0 votes {(full == 0).sum()}")

aucs, rho = [], []
for _ in range(1000):
    i = rng.choice(len(M), len(M))
    aucs.append(roc_auc_score(y, vote(M[i])))
    p = rng.permutation(len(M))
    rho.append(spearmanr(vote(M[p[:len(M) // 2]]), vote(M[p[len(M) // 2:]]))[0])
print(f"resampling MODELS (pairs fixed): AUC mean {np.mean(aucs):.3f}, sd {np.std(aucs):.3f}, 95% [{np.percentile(aucs, 2.5):.2f},{np.percentile(aucs, 97.5):.2f}]")
print(f"split-half Spearman of pair rankings (225 vs 225 models): {np.mean(rho):.3f}; Spearman-Brown for 450: {2 * np.mean(rho) / (1 + np.mean(rho)):.3f}")

# how many models are needed? AUC spread vs ensemble size
for k in (50, 150, 450):
    a = [roc_auc_score(y, vote(M[rng.choice(len(M), k)])) for _ in range(500)]
    print(f"  {k:3} models: AUC {np.mean(a):.3f} +/- {np.std(a):.3f}")
# per-model quality: does any single model predict well?
pm = [roc_auc_score(y, np.nan_to_num(M[i], nan=0)) for i in range(len(M))]
print(f"single models: AUC median {np.median(pm):.3f}, max {np.max(pm):.3f}; ensemble {roc_auc_score(y, full):.3f}")
