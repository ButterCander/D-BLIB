import json, glob, os
import numpy as np, pandas as pd
D = "mteapy/src/mteapy/data"
meta = pd.read_csv(f"{D}/task_metadata.tsv", sep="\t").set_index("ID")
tasks = json.load(open("out/task_genes.json"))
EXP = ["X138", "X140", "X141"]          # cholesterol, mevalonate, farnesyl-PP synthesis: fixed before looking
print("Positive control: statins in pancreatic cancer cells (GSE149566). Expected: cholesterol / mevalonate / farnesyl-PP synthesis tasks UP (SREBP2 feedback)")
rows = []
for f in sorted(glob.glob("out/tidex_statin_*.tsv")):
    name = os.path.basename(f)[len("tidex_statin_"):-4]
    r = pd.read_csv(f, sep="\t").set_index("task_id"); e = pd.read_csv(f"out/tideex_statin_{name}.tsv", sep="\t").set_index("task_id")
    s = r["score"].dropna(); n = len(s); rk = s.rank(ascending=False)
    # naive mean of masked logFC over task genes
    ln, tr = name.split("_", 1)
    de = pd.read_csv(f"out/statin_de_{name}.tsv", sep="\t")
    lf = de.set_index("gene")["logFC"].where(de.set_index("gene")["adj.P.Val"] < .05, 0.0)
    rows.append(dict(contrast=name, n_tasks=n, **{f"{t}_score": s.get(t, np.nan) for t in EXP}, **{f"{t}_rank": int(rk.get(t, np.nan)) if t in rk.index else -1 for t in EXP},
                     best_of_3_rank=int(min(rk.get(t, 999) for t in EXP)), n_up_p025=int(((r.score > 0) & (r.pvalue < .025)).sum()), n_down_p025=int(((r.score < 0) & (r.pvalue < .025)).sum()),
                     exp_p=",".join(f"{r.pvalue.get(t, np.nan):.3f}" for t in EXP)))
df = pd.DataFrame(rows).set_index("contrast")
pd.set_option("display.width", 250)
print(df[["n_tasks", "X138_score", "X140_score", "X141_score", "X138_rank", "X140_rank", "X141_rank", "best_of_3_rank", "exp_p", "n_up_p025", "n_down_p025"]].round(2).to_string())
# top tasks for one contrast
for name in ("MiaPaca2_Simvastatin_High", "PANC1_Simvastatin_High"):
    r = pd.read_csv(f"out/tidex_statin_{name}.tsv", sep="\t").set_index("task_id")
    r["desc"] = [meta.loc[i, "DESCRIPTION"][:44] for i in r.index]
    print(f"\n{name}: top 6 up / top 4 down by TIDE score")
    print(r.sort_values("score", ascending=False).head(6)[["desc", "score", "pvalue"]].round(3).to_string())
    print(r.sort_values("score").head(4)[["desc", "score", "pvalue"]].round(3).to_string())
# TIDE-essential
print("\nTIDE-essential: expected tasks (score, p)")
for f in sorted(glob.glob("out/tideex_statin_*High*.tsv")):
    name = os.path.basename(f)[len("tideex_statin_"):-4]; r = pd.read_csv(f, sep="\t").set_index("task_id")
    s = r["score"].dropna(); rk = s.rank(ascending=False)
    print(f"  {name:28} " + " ".join(f"{t}: {r.score.get(t, np.nan):+.2f} (rank {int(rk[t]) if t in rk.index else -1}/{len(s)}, p {r.pvalue.get(t, np.nan):.3f})" for t in EXP))
