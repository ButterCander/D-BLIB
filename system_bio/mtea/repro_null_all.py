"""Do the other four conditions give the same picture? authors' scores under zero-padded nulls of size N."""
import time, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
COND = {"PI103": "PI3Ki", "PD0325901": "MEKi", "5Z": "TAKi", "PI103+5Z": "PI3Ki-TAKi", "PI103+PD": "PI3Ki-MEKi"}
NS = [None, 4000, 5000, 6000, 8000, 16932]     # None = model genes with data only
for c, sheet in COND.items():
    r0 = pd.read_csv(f"out/de_{c}.tsv", sep="\t", index_col=0); r0["padj"] = r0["padj"].fillna(1)
    mg = r0.loc[r0.index.intersection(gm.index)].copy(); mg.index = [gm[i] for i in mg.index]
    mg["lfc"] = np.where(mg.padj < .05, mg.lfc_shrunk, 0.0)
    au = pd.read_excel("si/MOESM4.xlsx", sheet_name=sheet, header=None).iloc[2:, [0, 2, 3]]
    au.columns = ["tid", "sa", "pa"]; au["tid"] = "X" + au.tid.astype(str); au["pa"] = pd.to_numeric(au.pa); au["sa"] = pd.to_numeric(au.sa)
    line = []
    for N in NS:
        n = len(mg) if N is None else N
        e = pd.DataFrame({"lfc": mg.lfc.values}, index=mg.index)
        if n > len(e): e = pd.concat([e, pd.DataFrame({"lfc": np.zeros(n - len(e))}, index=[f"d{i}" for i in range(n - len(e))])])
        res = compute_TIDE(e, "lfc", ts, model, "absmax", 2000, 1, True, 1)
        nulls = {r.task_id: np.array(r.random_score_array.split(";"), dtype=float) for r in res.itertuples()}
        p = np.array([min((nulls[t] <= s).mean(), (nulls[t] >= s).mean()) if t in nulls else np.nan for t, s in zip(au.tid, au.sa)])
        ok = ~np.isnan(p); q = bh(np.clip(2 * p[ok], 1e-12, 1))
        lp = lambda x: -np.log10(np.clip(x, 1e-4, 1)); rm = np.sqrt(np.mean((lp(au.pa.values[ok]) - lp(p[ok])) ** 2))
        line.append(f"N={'model' if N is None else N}: {int((p[ok]<.025).sum())}/{int(q.__lt__(.05).sum())} (RMSE {rm:.2f})")
    qa = bh(np.clip(2 * au.pa.values, 1e-12, 1))
    print(f"{c:10} authors p<.025 {(au.pa<.025).sum()} / BH {(qa<.05).sum()} | nonzero model genes {(mg.lfc!=0).sum()} | " + " | ".join(line), flush=True)
print("done")
