"""Analyse syn_oof.npz: tuned CV table, paired cluster-bootstrap deltas (resample drug pairs), per-cell-line, figures."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from syn_final import *

oof = np.load("syn_oof.npz")
n, yb = len(y), y > 0
G = {g: np.where(pair.values == g)[0] for g in pair.unique()}
mean_pred = lambda sc, fs, mn: oof[key(sc, fs, mn)].mean(0)
auc = lambda yt, p: roc_auc_score(yt > 0, p)


def table(sc):
    print(f"== {sc}: mean +/- sd over {R} repeats of 5-fold, nested-tuned   (Pearson r | AUC)")
    print(f"{'features':22}" + "".join(f"{m:>32}" for m in GRIDS))
    for fs in FS:
        row = ""
        for mn in GRIDS:
            P = oof[key(sc, fs, mn)]
            r = [pearsonr(y, p)[0] for p in P]
            a = [auc(y, p) for p in P]
            row += f"  {np.mean(r):.3f}+/-{np.std(r):.3f} | {np.mean(a):.3f}+/-{np.std(a):.3f}"
        print(f"{fs:22}{row}")
    print()


def delta(pa, pb, B=2000):
    """pa - pb for Pearson r and AUC; 95% CI by resampling drug pairs (clusters)."""
    f = lambda s: (pearsonr(y[s], pa[s])[0] - pearsonr(y[s], pb[s])[0], auc(y[s], pa[s]) - auc(y[s], pb[s]))
    rng, ks, bs = np.random.default_rng(0), list(G), []
    for _ in range(B):
        bs.append(f(np.concatenate([G[k] for k in rng.choice(ks, len(ks))])))
    return f(np.arange(n)), np.percentile(bs, [2.5, 97.5], axis=0)


mp = lambda *a: mean_pred("pair", *a)
COMP = {
    "SVR drug+cell  vs  RF Boolean (paper)": (mp("drug+cell", "SVR"), mp("Boolean", "RF")),
    "Ridge drug+cell  vs  RF Boolean (paper)": (mp("drug+cell", "Ridge"), mp("Boolean", "RF")),
    "SVR target+cell  vs  RF Boolean (paper)": (mp("target+cell", "SVR"), mp("Boolean", "RF")),
    "SVR Boolean+drug+cell  vs  SVR drug+cell": (mp("Boolean+drug+cell", "SVR"), mp("drug+cell", "SVR")),
    "Ridge Boolean+drug+cell  vs  Ridge drug+cell": (mp("Boolean+drug+cell", "Ridge"), mp("drug+cell", "Ridge")),
    "Ensemble(SVR drug+cell, RF Boolean)  vs  SVR drug+cell": ((mp("drug+cell", "SVR") + mp("Boolean", "RF")) / 2, mp("drug+cell", "SVR")),
}

if __name__ == "__main__":
    table("random")
    table("pair")
    print("== paired differences, pair-grouped CV, repeat-averaged predictions; [95% CI, resampling drug pairs]")
    res = {}
    for name, (a, b) in COMP.items():
        d, ci = delta(a, b)
        res[name] = (d, ci)
        print(f"{name:58} dr={d[0]:+.3f} [{ci[0][0]:+.3f},{ci[1][0]:+.3f}]   dAUC={d[1]:+.3f} [{ci[0][1]:+.3f},{ci[1][1]:+.3f}]")
    print("\n== per cell line, pair-grouped CV (Pearson r | AUC)")
    print(f"{'cell line':12}{'RF Boolean (paper)':>22}{'SVR drug+cell':>22}{'SVR Boolean+drug+cell':>26}")
    for c in sorted(cell.unique()):
        s = (cell == c).values
        cols = [mp("Boolean", "RF"), mp("drug+cell", "SVR"), mp("Boolean+drug+cell", "SVR")]
        print(f"{c:12}" + "".join(f"{pearsonr(y[s], p[s])[0]:>14.3f} | {auc(y[s], p[s]):.3f}" for p in cols))

    os.makedirs("fig", exist_ok=True)
    col = {"RF": "#9aa0a6", "Ridge": "#4c78a8", "SVR": "#e45756"}
    fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
    for k, (metric, fn) in enumerate({"Pearson r": pearsonr, "AUC": auc}.items()):
        for j, mn in enumerate(GRIDS):
            v = [[(fn(y, p)[0] if k == 0 else fn(y, p)) for p in oof[key("pair", fs, mn)]] for fs in FS]
            ax[k].bar(np.arange(len(FS)) + (j - 1) * .27, np.mean(v, 1), .27, yerr=np.std(v, 1), color=col[mn], label=mn, capsize=2)
        ax[k].set_xticks(range(len(FS)), [f.replace("+", "+\n") for f in FS], fontsize=8)
        ax[k].set_title(metric + " (new drug pairs, nested CV)"), ax[k].spines[["top", "right"]].set_visible(False)
    ax[0].axhline(.403, ls="--", c="k", lw=1), ax[0].text(-.45, .545, "dashed = paper baseline 0.403", fontsize=8)
    ax[1].axhline(.674, ls="--", c="k", lw=1), ax[1].text(-.45, .795, "dashed = paper baseline 0.674", fontsize=8), ax[1].legend(frameon=False, fontsize=8, loc="upper right")
    ax[0].set_ylim(0), ax[1].set_ylim(.5)
    plt.tight_layout(), plt.savefig("fig/fig1_models.png", dpi=200), plt.close()

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.2), sharey=True)
    for k, metric in enumerate(["Pearson r", "AUC"]):
        for i, (name, (d, ci)) in enumerate(res.items()):
            ax[k].errorbar(d[k], i, xerr=[[d[k] - ci[0][k]], [ci[1][k] - d[k]]], fmt="o", c="#e45756" if ci[0][k] > 0 else "#555", capsize=3)
        ax[k].axvline(0, c="k", lw=1), ax[k].set_title("difference in " + metric + " (95% CI)")
        ax[k].spines[["top", "right"]].set_visible(False)
    ax[0].set_yticks(range(len(res)), [k.replace("  vs  ", "\n  vs  ") for k in res], fontsize=7), ax[0].invert_yaxis()
    plt.tight_layout(), plt.savefig("fig/fig2_deltas.png", dpi=200), plt.close()
    print("\nfigures: fig/fig1_models.png fig/fig2_deltas.png")
