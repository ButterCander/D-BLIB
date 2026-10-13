"""Slide figures for the Flobak 2023 presentation (1664 px wide at 200 dpi = full content width of a 1920x1080 slide)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import rankdata
from ags_eval import *

BG, BLUE, ORANGE, INK, GRAY, MUTED = "#F7F6F2", "#0072B2", "#D55E00", "#13243A", "#8A94A3", "#4A5568"
plt.rcParams.update({"font.family": "Arial", "font.size": 13, "axes.edgecolor": GRAY, "axes.linewidth": 1.2,
                     "axes.spines.top": False, "axes.spines.right": False, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK, "text.color": INK, "figure.facecolor": BG, "axes.facecolor": BG})
W, DPI = 1664 / 200, 200
pos, neg = np.where(y)[0], np.where(~y)[0]


def boot(s, B=2000, seed=0):
    rng, out = np.random.default_rng(seed), []
    for _ in range(B):
        i = np.concatenate([rng.choice(pos, len(pos)), rng.choice(neg, len(neg))])
        out.append((roc_auc_score(y[i], s[i]), average_precision_score(y[i], s[i])))
    return np.percentile(out, [2.5, 97.5], axis=0)


def stat(s):
    lo, hi = boot(s)
    return roc_auc_score(y, s), average_precision_score(y, s), lo, hi


FOLD = {"Parameters only": "link-only", "Topology only": "topology-only", "Topology + parameters": "topo-and-link"}
exn = lambda f, k: load(f, k, 150, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
# paper rule: calibrated + beta * random-model excess with beta = -1, negative = synergy  ->  score = random - calibrated (higher = synergy)
EX = {k: exn(f, "rand") - exn(f, "ss") for k, f in FOLD.items()}
V = {f: frac(load(f, "ss", 150, level="modelwise")).fillna(0).values for f in FOLD.values()}

# ---- Fig A: paper-style AUC with bootstrap CI --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(W, 3.1))
print("FIG A")
for r, (k, s) in enumerate(EX.items()):
    a, ap, lo, hi = stat(s)
    yy = 2 - r
    ax.errorbar(a, yy, xerr=[[a - lo[0]], [hi[0] - a]], fmt="o", color=BLUE, ms=13, capsize=7, lw=3.5, capthick=3)
    ax.text(1.03, yy, f"{a:.2f}   [{lo[0]:.2f}, {hi[0]:.2f}]", va="center", fontsize=14, fontweight="bold")
    print(f"  {k:24} AUC {a:.3f} [{lo[0]:.2f},{hi[0]:.2f}]  PR {ap:.3f} [{lo[1]:.2f},{hi[1]:.2f}]")
ax.axvline(.5, ls="--", color=GRAY, lw=2), ax.text(.51, 2.62, "chance", ha="left", color=MUTED, fontsize=13)
ax.set_yticks([2, 1, 0], list(EX)), ax.set_xlim(.2, 1.32), ax.set_ylim(-.6, 2.9)
ax.set_xticks([.2, .4, .6, .8, 1.0]), ax.set_xlabel("ROC AUC (dot) with 95% bootstrap interval")
ax.tick_params(axis="y", length=0, labelsize=14), ax.spines["left"].set_visible(False)
plt.tight_layout(), plt.savefig("fig/ags_A_auc_ci.png", dpi=DPI), plt.close()

# ---- Fig B: single models vs ensemble ----------------------------------------------------------------------
M = pd.read_csv(f"{A}/results/link-only/cascade_2.0_ss_150sim_fixpoints_bliss_model_predictions.tab", sep="\t", index_col=0)
M.columns = [pk(c) for c in M.columns]
M = M[pairs].values
vote = lambda m: np.nansum(m, 0) / np.maximum((~np.isnan(m)).sum(0), 1)
single = np.array([roc_auc_score(y, np.nan_to_num(M[i], nan=0)) for i in range(len(M))])
full = roc_auc_score(y, vote(M))
rng = np.random.default_rng(0)
ks = [10, 25, 50, 100, 150, 250, 350]
curve = [[roc_auc_score(y, vote(M[rng.choice(len(M), k, replace=False)])) for _ in range(400)] for k in ks]
print(f"FIG B  single-model AUC median {np.median(single):.3f} (min {single.min():.2f}, max {single.max():.2f}); ensemble {full:.3f}; "
      f"missing {np.isnan(M).mean():.1%}")
for k, c in zip(ks, curve):
    print(f"  k={k:3}: {np.mean(c):.3f} +/- {np.std(c):.3f}")
fig, (a1, a2) = plt.subplots(1, 2, figsize=(W, 3.3), gridspec_kw={"width_ratios": [1, 1.05]})
a1.hist(single, bins=np.linspace(.3, .9, 25), color=GRAY, edgecolor=BG)
a1.axvline(full, color=ORANGE, lw=4), a1.text(full - .01, a1.get_ylim()[1] * .86, f"ensemble\n{full:.2f}", color=ORANGE, ha="right", fontsize=14, fontweight="bold")
a1.axvline(.5, ls="--", color=INK, lw=2), a1.text(.44, a1.get_ylim()[1] * .55, f"single models\nmedian {np.median(single):.2f}", ha="right", fontsize=13)
a1.set_xlabel("AUC of one model"), a1.set_ylabel("Number of models"), a1.set_xlim(.2, .9)
m, s = np.array([np.mean(c) for c in curve]), np.array([np.std(c) for c in curve])
a2.fill_between(ks, m - s, m + s, color=BLUE, alpha=.22, lw=0), a2.plot(ks, m, "-o", color=BLUE, lw=3.5, ms=8)
a2.plot([len(M)], [full], "o", color=ORANGE, ms=14), a2.text(len(M) - 12, full - .075, "all 450", ha="right", color=ORANGE, fontsize=14, fontweight="bold")
a2.axhline(.5, ls="--", color=GRAY, lw=2), a2.text(465, .515, "chance", ha="right", color=MUTED, fontsize=13)
a2.set_ylim(.4, .92), a2.set_xlim(0, 470), a2.set_xlabel("Models in the ensemble"), a2.set_ylabel("AUC (mean ± SD)")
plt.tight_layout(), plt.savefig("fig/ags_B_ensemble.png", dpi=DPI), plt.close()

# ---- Fig C: consensus scores, AUC and PR AUC with CI ---------------------------------------------------------
vs, rl = list(V.values()), list(EX.values())
rows = {"Paper rule, parameters only": rl[0], "Paper rule, topology + parameters": rl[2], "Vote fraction, parameters only": vs[0],
        "Consensus of 6 scores": np.mean([rankdata(v) for v in (*rl, *vs)], 0)}
pi = np.array(["PI" in p for p in pairs], float)
print("FIG C")
fig, axs = plt.subplots(1, 2, figsize=(W, 3.4), sharey=True)
for r, (name, s) in enumerate(rows.items()):
    a, ap, lo, hi = stat(s)
    yy = len(rows) - 1 - r
    col = ORANGE if r == len(rows) - 1 else BLUE
    for ax, val, l, h in ((axs[0], a, lo[0], hi[0]), (axs[1], ap, lo[1], hi[1])):
        ax.errorbar(val, yy, xerr=[[val - l], [h - val]], fmt="o", color=col, ms=11, capsize=6, lw=3, capthick=2.5)
        ax.text(h + .025, yy, f"{val:.2f}", va="center", fontsize=13, fontweight="bold", color=col)
    print(f"  {name:38} AUC {a:.3f} [{lo[0]:.2f},{hi[0]:.2f}]  PR AUC {ap:.3f} [{lo[1]:.2f},{hi[1]:.2f}]")
pa, pp = roc_auc_score(y, pi), average_precision_score(y, pi)
print(f"  trivial 'pair contains PI': AUC {pa:.3f}, PR AUC {pp:.3f}; chance PR {y.mean():.3f}")
for ax, base, chance, t in ((axs[0], pa, .5, "ROC AUC"), (axs[1], pp, y.mean(), "PR AUC")):
    ax.axvline(chance, ls="--", color=GRAY, lw=2), ax.axvline(base, ls=":", color=INK, lw=3)
    ax.set_title(t, loc="left", fontsize=15, fontweight="bold"), ax.set_ylim(-.6, len(rows) - .3)
axs[0].set_xlim(.3, 1.08), axs[1].set_xlim(0, 1.08)
axs[0].set_yticks(range(len(rows)), list(rows)[::-1], fontsize=13), axs[0].tick_params(axis="y", length=0)
axs[1].tick_params(left=False, labelleft=False)
for ax in axs:
    ax.spines["left"].set_visible(False), ax.set_xlabel("estimate with 95% interval", fontsize=12)
from matplotlib.lines import Line2D
fig.legend([Line2D([], [], ls="--", color=GRAY, lw=2), Line2D([], [], ls=":", color=INK, lw=3)],
           ["chance", "trivial rule: pair contains PI"], loc="upper center", ncol=2, frameon=False, fontsize=13)
plt.tight_layout(rect=(0, 0, 1, .92)), plt.savefig("fig/ags_C_consensus.png", dpi=DPI), plt.close()
