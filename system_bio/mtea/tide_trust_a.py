"""Trust check A: does the GPR/reaction projection of TIDE add information over a naive mean of the (masked) logFC of the task's genes?"""
import json
import numpy as np, pandas as pd
from scipy import stats
D = "mteapy/src/mteapy/data"
tasks = json.load(open("out/task_genes.json"))
meta = pd.read_csv(f"{D}/task_metadata.tsv", sep="\t").set_index("ID")
CONDS = ["PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]
print("[A] GSE285616, 24 h: TIDE score (model-projected) vs naive mean masked logFC of the task's genes, across 187 tasks")
rows = []
for c in CONDS:
    de = pd.read_csv(f"out/de_{c}.tsv", sep="\t", index_col=0)
    lf = de["lfc_shrunk"].where(de["padj"].fillna(1) < .05, 0.0)
    naive = pd.Series({t: lf.reindex(g).dropna().mean() for t, g in tasks.items() if lf.reindex(g).dropna().size >= 1})
    tide = pd.read_csv(f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")["score"]
    j = pd.concat([tide, naive], axis=1, keys=["tide", "naive"]).dropna()
    top = 15
    a = set(j.tide.abs().nlargest(top).index); b = set(j.naive.abs().nlargest(top).index)
    rows.append((c, len(j), stats.pearsonr(j.tide, j.naive)[0], stats.spearmanr(j.tide, j.naive)[0], len(a & b)))
    print(f"    {c:10} n={len(j)} Pearson {rows[-1][2]:.2f} Spearman {rows[-1][3]:.2f} | top-{top} |score| overlap {rows[-1][4]}/{top}")
print("    tasks where the sign of TIDE and the naive mean disagree (|TIDE|>0.2 and |naive|>0.2):")
for c in CONDS[:2]:
    de = pd.read_csv(f"out/de_{c}.tsv", sep="\t", index_col=0)
    lf = de["lfc_shrunk"].where(de["padj"].fillna(1) < .05, 0.0)
    naive = pd.Series({t: lf.reindex(g).dropna().mean() for t, g in tasks.items()})
    tide = pd.read_csv(f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")["score"]
    j = pd.concat([tide, naive], axis=1, keys=["tide", "naive"]).dropna()
    bad = j[(j.tide.abs() > .2) & (j.naive.abs() > .2) & (np.sign(j.tide) != np.sign(j.naive))]
    print(f"    {c}: {len(bad)} of {int(((j.tide.abs()>.2)&(j.naive.abs()>.2)).sum())}")
