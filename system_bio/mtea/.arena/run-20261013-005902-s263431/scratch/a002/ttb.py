"""Trust checks B-E + PDAC application.
B: task-level replication of TIDE between two independent AGS experiments (GSE285616 half dose vs Tsirvouli full dose), matched vs unmatched drugs.
C: replication across four PDAC lines and across time (trametinib).  D: cross-cancer agreement (PDAC trametinib vs AGS MEKi, vs AGS PI3Ki as a specificity control).
E: consensus PDAC tasks and their DepMap pancreatic dependency."""
import json, glob
import numpy as np, pandas as pd
from scipy import stats
D = "mteapy/src/mteapy/data"
meta = pd.read_csv(f"{D}/task_metadata.tsv", sep="\t").set_index("ID")
ess = pd.read_csv(f"{D}/HumanGEM_essential_genes_matrix.tsv", sep="\t", index_col=0).astype(bool)
hg = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneID").set_index("geneID")["geneSymbol"]
def load(path):
    r = pd.read_csv(path, sep="\t"); return r.set_index("task_id")["score"]
GS = {"PI": "PI103", "PD": "PD0325901", "OXO": "5Z", "PIOXO": "PI103+5Z", "PIPD": "PI103+PD"}
T = {k: load(f"out/tidex_tsir_{k}_24h.tsv") for k in GS}
G = {k: load(f"out/tide_p10000_{v}.tsv") for k, v in GS.items()}
TE = {k: load(f"out/tideex_tsir_{k}_24h.tsv") for k in GS}
GE = {k: load(f"out/tidee_{v}.tsv") for k, v in GS.items()}
def rho(a, b, m=0.0):
    j = pd.concat([a, b], axis=1).dropna(); return stats.spearmanr(j.iloc[:, 0], j.iloc[:, 1])[0], len(j)

print("[B] task-level replication, independent AGS experiments: Spearman of TIDE scores, rows = GSE285616 (half dose, 24 h), cols = Tsirvouli (full dose, 24 h)")
for lab, X, Y in (("TIDE", G, T), ("TIDE-essential", GE, TE)):
    M = pd.DataFrame({c: {r: rho(X[r], Y[c])[0] for r in GS} for c in GS})
    print(f"  {lab}\n" + M.round(2).to_string().replace("\n", "\n  "))
    d = np.diag(M.values); off = M.values[~np.eye(5, dtype=bool)]
    print(f"    matched drug: mean {d.mean():.2f} (min {d.min():.2f}); unmatched: mean {off.mean():.2f} (max {off.max():.2f})")

pl = ["CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"]
P = {(l, k): load(f"out/tidex_pdac_{l}_{k}.tsv") for l in pl for k in ("h3", "h12", "h24")}
PE = {(l, k): load(f"out/tideex_pdac_{l}_{k}.tsv") for l in pl for k in ("h3", "h12", "h24")}
print("\n[C] PDAC trametinib: Spearman of TIDE scores between cell lines")
for k in ("h12", "h24"):
    M = pd.DataFrame({b: {a: rho(P[(a, k)], P[(b, k)])[0] for a in pl} for b in pl})
    off = M.values[~np.eye(4, dtype=bool)]
    print(f"  {k}: mean off-diagonal {off.mean():.2f} (range {off.min():.2f}..{off.max():.2f})\n  " + M.round(2).to_string().replace("\n", "\n  "))
print("  within line, 12 h vs 24 h:", {l: round(rho(P[(l, 'h12')], P[(l, 'h24')])[0], 2) for l in pl})
print("  TIDE-essential, 24 h between lines: mean off-diag", round(np.mean([rho(PE[(a, 'h24')], PE[(b, 'h24')])[0] for i, a in enumerate(pl) for b in pl[i+1:]]), 2))

print("\n[D] cross-cancer: PDAC trametinib (24 h) vs AGS single drugs (Spearman of TIDE scores)")
rows = []
for l in pl:
    rows.append([l] + [round(rho(P[(l, "h24")], G[k])[0], 2) for k in ("PD", "PI", "OXO")] + [round(rho(P[(l, "h24")], T[k])[0], 2) for k in ("PD", "PI", "OXO")])
print(pd.DataFrame(rows, columns=["line", "GSE285616 MEKi", "PI3Ki", "TAKi", "Tsirvouli MEKi", "PI3Ki", "TAKi"]).to_string(index=False))

print("\n[E] consensus PDAC tasks at 24 h: |score|>=0.15 with the same sign in >=3 of 4 lines")
S = pd.DataFrame({l: P[(l, "h24")] for l in pl}).dropna(thresh=3)
cons = []
for t, r in S.iterrows():
    v = r.dropna(); s = np.sign(v[v.abs() >= .15])
    for sg in (-1, 1):
        if (s == sg).sum() >= 3: cons.append((t, sg, v.mean()))
cdf = pd.DataFrame(cons, columns=["task", "sign", "mean"]).set_index("task")
ags = pd.concat([G["PD"], T["PD"]], axis=1, keys=["ags_g", "ags_t"])
cdf = cdf.join(ags); cdf["desc"] = [meta.loc[i, "DESCRIPTION"][:46] if i in meta.index else "" for i in cdf.index]
cdf["agree_ags"] = (np.sign(cdf["mean"]) == np.sign(cdf.ags_g)) & (np.sign(cdf["mean"]) == np.sign(cdf.ags_t))
print(f"  {len(cdf)} consensus tasks ({int((cdf.sign<0).sum())} down, {int((cdf.sign>0).sum())} up); same direction in both AGS MEKi datasets: {int(cdf.agree_ags.sum())}")

try:
    cr = pd.read_csv("external/crispr_pancreas_rows.tsv", sep="\t", index_col=0)
    cr = cr.loc[:, ~cr.columns.duplicated()]
    def dep(t):
        g = [hg[x] for x in ess.index[ess[t]] if x in hg.index]; v = cr.reindex(columns=g).dropna(axis=1, how="all")
        if v.shape[1] == 0: return np.nan, np.nan, 0
        mn = v.min(axis=1)                      # per line: strongest dependency among the task's essential genes
        return mn.median(), (mn <= -0.5).mean(), v.shape[1]
    out = [(t,) + dep(t) for t in cdf.index if t in ess.columns]
    dd = pd.DataFrame(out, columns=["task", "min_chronos_median", "frac_lines_dep", "n_genes"]).set_index("task")
    cdf = cdf.join(dd)
    allt = pd.DataFrame([(t,) + dep(t) for t in ess.columns], columns=["task", "m", "f", "n"]).set_index("task")
    print(f"  median min-Chronos in pancreatic lines: consensus tasks {cdf.min_chronos_median.median():+.2f} vs all tasks {allt.m.median():+.2f}")
    pd.set_option("display.width", 200)
    dn = cdf[cdf["sign"] < 0].sort_values("mean").head(14)
    up = cdf[cdf["sign"] > 0].sort_values("mean", ascending=False).head(8)
    cols = ["desc", "mean", "ags_g", "ags_t", "min_chronos_median", "frac_lines_dep", "n_genes"]
    print("  DOWN (top 14 by mean score):\n" + dn[cols].round(2).to_string().replace("\n", "\n  "))
    print("  UP (top 8):\n" + up[cols].round(2).to_string().replace("\n", "\n  "))
    cdf.to_csv("out/pdac_consensus_tasks.tsv", sep="\t")
except FileNotFoundError:
    print("  DepMap pancreas rows not found")
