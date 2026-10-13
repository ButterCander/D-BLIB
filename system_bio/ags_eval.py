"""Flobak 2023 (CASCADE 2.0, AGS): reproduce ROC/PR AUC of drug-synergy predictions from the published result tables."""
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

A = r"D:\Cander\School\3-2\system_bio\ags"
pk = lambda s: frozenset(s.replace("[", "").replace("]", "").split("-"))
obs = {frozenset(l.strip().split("~")) for l in open(f"{A}/data/observed_synergies_cascade_2.0") if l.strip()}


def load(folder, kind, nsim, metric="bliss", level="ensemblewise"):
    d = pd.read_csv(f"{A}/results/{folder}/cascade_2.0_{kind}_{nsim}sim_fixpoints_{metric}_{level}_synergies.tab", sep="\t", index_col=0)
    d.index = [pk(i) for i in d.index]
    return d


pairs = list(load("link-only", "ss", 150, level="modelwise").index)
y = np.array([p in obs for p in pairs])
assert y.sum() == 6 and len(pairs) == 153


def frac(d):  # fraction of models predicting synergy
    return (d.Synergies / (d.Synergies + d["Non-synergies"])).reindex(pairs)


def scores(folder, metric="bliss", nsim=150):
    ens_ss, ens_r = (load(folder, k, nsim, metric)["Response excess over subset"].reindex(pairs) for k in ("ss", "rand"))
    p_ss, p_r = (frac(load(folder, k, nsim, metric, "modelwise")) for k in ("ss", "rand"))
    e = 1e-3
    return {"ensemble excess (ss)": ens_ss, "ensemble excess ss-rand": ens_ss - ens_r,
            "model vote fraction (ss)": p_ss, "vote fraction ss-rand": p_ss - p_r,
            "vote log2 fold-change ss/rand": np.log2((p_ss + e) / (p_r + e))}


def ci(s, B=5000, seed=0):
    """stratified bootstrap (resample positives and negatives separately) and permutation p-value for AUC."""
    rng = np.random.default_rng(seed)
    pos, neg = np.where(y)[0], np.where(~y)[0]
    boot = []
    for _ in range(B):
        i = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        boot.append((roc_auc_score(y[i], s[i]), average_precision_score(y[i], s[i])))
    perm = np.mean([roc_auc_score(rng.permutation(y), s) >= roc_auc_score(y, s) for _ in range(B)])
    return np.percentile(boot, [2.5, 97.5], axis=0), perm


if __name__ == "__main__":
    print(f"153 pairs, 6 observed synergies; random-model AUC = 0.50, PR AUC baseline = {y.mean():.3f}\n")
    for folder in ("link-only", "topology-only", "topo-and-link"):
        for metric in ("bliss", "hsa"):
            print(f"== {folder}, {metric}, 150 evolution runs")
            for name, s in scores(folder, metric).items():
                v = s.fillna(s.min()).values
                (lo, hi), p = ci(v, 2000)
                print(f"  {name:32} AUC {roc_auc_score(y, v):.3f} [{lo[0]:.2f},{hi[0]:.2f}]  PR AUC {average_precision_score(y, v):.3f} [{lo[1]:.2f},{hi[1]:.2f}]  perm p={p:.3f}")
