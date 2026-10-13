"""Which null universe reproduces the authors' PI3Ki p-values from the authors' own scores?
Zero-padding the gene vector (genes tested but not significant, or untested genes) narrows the null for the gene-label permutation."""
import time, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
r0 = pd.read_csv("out/de_PI103.tsv", sep="\t", index_col=0); r0["padj"] = r0["padj"].fillna(1)
cnt_rows = pd.read_csv("GSE285616_counts.csv.gz", usecols=["Symbol"]).Symbol.nunique()
print(f"genes in count file: {cnt_rows}; DE-tested: {len(r0)}", flush=True)
mg = r0.loc[r0.index.intersection(gm.index)].copy(); mg.index = [gm[i] for i in mg.index]
mg["lfc"] = np.where(mg.padj < .05, mg.lfc_shrunk, 0.0)
au = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki", header=None).iloc[2:, [0, 2, 3]]
au.columns = ["tid", "sa", "pa"]; au["tid"] = "X" + au.tid.astype(str); au["pa"] = pd.to_numeric(au.pa); au["sa"] = pd.to_numeric(au.sa)
n_nonmodel_tested = len(r0) - len(mg)
print(f"model genes tested {len(mg)} (nonzero {(mg.lfc!=0).sum()}), non-model tested {n_nonmodel_tested}", flush=True)
UNIV = {f"N={n}": n for n in (2443, 5000, 16932)}
for name, N in UNIV.items():
    e = pd.DataFrame({"lfc": mg.lfc.values}, index=mg.index)
    extra = max(0, N - len(e))
    if extra:
        e = pd.concat([e, pd.DataFrame({"lfc": np.zeros(extra)}, index=[f"dummy{i}" for i in range(extra)])])
    t0 = time.time()
    res = compute_TIDE(e, "lfc", ts, model, "absmax", 1000, 1, True, 7)
    nulls = {r.task_id: np.array(r.random_score_array.split(";"), dtype=float) for r in res.itertuples()}
    p = []
    for r in au.itertuples():
        if r.tid in nulls:
            n = nulls[r.tid]; p.append(min((n <= r.sa).mean(), (n >= r.sa).mean()))
        else: p.append(np.nan)
    au["p_u"] = p
    ok = au.p_u.notna()
    q = bh(np.clip(2 * au.p_u[ok].values, 1e-12, 1))
    rho = stats.spearmanr(au.pa[ok], au.p_u[ok])[0]
    lp = lambda x: -np.log10(np.clip(x, 1e-4, 1)); rm = float(np.sqrt(np.mean((lp(au.pa[ok].values) - lp(au.p_u[ok].values)) ** 2)))
    print(f"{name:34} N={N:6d} nonzero frac {(e.lfc!=0).mean():.3f} | authors' scores under this null: p<.025 {(au.p_u[ok]<.025).sum():3d} (authors' own p: {(au.pa<.025).sum()}) BH {(q<.05).sum():3d} (authors 59) | Spearman(p_authors, p_here) {rho:.2f} | RMSE(-log10 p) {rm:.2f} | {time.time()-t0:.0f}s", flush=True)
print("done")
