import pandas as pd, numpy as np, glob, os
P = r"D:/Cander/School/3-2/system_bio/mtea/out/"
f = pd.read_csv(P + "tidex_statin_PANC1_Simvastatin_High.tsv", sep="\t")
print(f.columns.tolist(), f.shape)
desc = {}
for cell in ("PANC1", "MiaPaca2"):
    for st in ("Fluvastatin", "Lovastatin", "Simvastatin"):
        for dose in ("High", "Low"):
            d = pd.read_csv(P + f"tidex_statin_{cell}_{st}_{dose}.tsv", sep="\t").set_index("task_id")
            s = d["score"].dropna()
            r_pos = s.rank(ascending=False)
            r_abs = s.abs().rank(ascending=False)
            out = []
            for t in ("X138", "X140", "X141"):
                if t in s.index:
                    out.append(f"{t}: score {s[t]:+.2f} rank_pos {int(r_pos[t])} rank_abs {int(r_abs[t])} p {d.loc[t,'pvalue']:.4f}")
                else:
                    out.append(f"{t}: NA")
            print(cell, st, dose, "n tasks", len(s), "n nonzero", int((s != 0).sum()), " | ", " | ".join(out))
# which tasks are these
t = pd.read_csv(P + "tidex_statin_PANC1_Simvastatin_High.tsv", sep="\t").set_index("task_id")
print(t.loc[["X138", "X140", "X141"]].to_string())
print(t.sort_values("score", ascending=False).head(8)[["score", "pvalue"] + [c for c in t.columns if c not in ("score", "pvalue")][:2]].to_string())
