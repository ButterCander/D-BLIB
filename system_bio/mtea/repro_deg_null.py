"""Null built from the AUTHORS' own DEG lists (their actual LFC values, |LFC|>=1 & padj<.05, all other entries 0) zero-padded to N;
authors' scores from Supplementary Data 3 evaluated against it."""
import numpy as np, pandas as pd
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
gids = sorted(set(pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").geneID))
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
print(f"Human-GEM genes: {len(gids)}", flush=True)
for sh in ["PI3Ki", "MEKi", "TAKi", "PI3Ki-TAKi", "PI3Ki-MEKi"]:
    deg = pd.read_excel("si/MOESM2.xlsx", sheet_name=sh)
    lf = deg[deg.geneID.isin(gids)].drop_duplicates("geneID").set_index("geneID")["log2FoldChange"]
    a = pd.read_excel("si/MOESM4.xlsx", sheet_name=sh, header=None).iloc[2:, [0, 2, 3]]; a.columns = ["id", "s", "p"]
    a["tid"] = "X" + a.id.astype(str); a["s"] = pd.to_numeric(a.s); a["p"] = pd.to_numeric(a.p)
    line = []
    for N in (len(gids), 5000, 8000, 12000, 20000):
        e = pd.DataFrame({"lfc": 0.0}, index=gids); e.loc[lf.index, "lfc"] = lf.values
        if N > len(e): e = pd.concat([e, pd.DataFrame({"lfc": np.zeros(N - len(e))}, index=[f"d{i}" for i in range(N - len(e))])])
        res = compute_TIDE(e, "lfc", ts, model, "absmax", 2000, 1, True, 1)
        nulls = {r.task_id: np.array(r.random_score_array.split(";"), dtype=float) for r in res.itertuples()}
        p = np.array([min((nulls[t] <= s).mean(), (nulls[t] >= s).mean()) for t, s in zip(a.tid, a.s)])
        q = bh(np.clip(2 * p, 1e-12, 1)); lp = lambda v: -np.log10(np.clip(v, 1e-4, 1))
        line.append(f"N={N}: {int((p<.025).sum())}/{int((q<.05).sum())} (RMSE {np.sqrt(np.mean((lp(a.p.values)-lp(p))**2)):.2f})")
    qa = bh(np.clip(2 * a.p.values, 1e-12, 1))
    print(f"{sh:11} model-gene DEGs {len(lf):3d} | authors {(a.p<.025).sum()}/{(qa<.05).sum()} | " + " | ".join(line), flush=True)
print("done")
