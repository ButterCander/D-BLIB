import sys, time, numpy as np, pandas as pd
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
ROOT = r"D:/Cander/School/3-2/system_bio/mtea/"
D = ROOT + "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
cond, sheet = sys.argv[1], sys.argv[2]
r0 = pd.read_csv(ROOT + f"out/de_{cond}.tsv", sep="\t", index_col=0); r0["padj"] = r0["padj"].fillna(1)
mg = r0.loc[r0.index.intersection(gm.index)].copy(); mg.index = [gm[i] for i in mg.index]
mg = mg[~mg.index.duplicated()]
mg["lfc"] = np.where(mg.padj < .05, mg.lfc_shrunk, 0.0)
raw = pd.read_excel(ROOT + "si/MOESM4.xlsx", sheet_name=sheet, header=None).iloc[2:, :4]; raw.columns = ["tid","desc","sa","pa"]
raw = raw.dropna(subset=["tid"]); raw["tid"] = "X" + raw.tid.astype(float).astype(int).astype(str)
raw["pa"] = pd.to_numeric(raw.pa); raw["sa"] = pd.to_numeric(raw.sa)
print(f"{cond}/{sheet}: model genes with data {len(mg)}, nonzero {(mg.lfc!=0).sum()} ({(mg.lfc!=0).mean():.2f}); authors p<.025 {(raw.pa<.025).sum()}, BH {(bh(np.clip(2*raw.pa.values,1e-12,1))<.05).sum()}", flush=True)
for N in [int(x) for x in sys.argv[3].split(",")]:
    e = pd.DataFrame({"lfc": mg.lfc.values}, index=mg.index)
    extra = max(0, N - len(e))
    if extra: e = pd.concat([e, pd.DataFrame({"lfc": np.zeros(extra)}, index=[f"dummy{i}" for i in range(extra)])])
    t0 = time.time()
    res = compute_TIDE(e, "lfc", ts, model, "absmax", int(sys.argv[4]), 1, True, 7)
    nulls = {r.task_id: np.array(r.random_score_array.split(";"), dtype=float) for r in res.itertuples()}
    p = []
    for r in raw.itertuples():
        n = nulls.get(r.tid)
        p.append(np.nan if n is None else min((n <= r.sa).mean(), (n >= r.sa).mean()))
    raw["pu"] = p; ok = raw.pu.notna()
    q = bh(np.clip(2*raw.pu[ok].values, 1e-12, 1))
    print(f"N={N:6d} nonzero frac {(e.lfc!=0).mean():.3f} | authors' scores under this null: p<.025 {(raw.pu[ok]<.025).sum()}, BH<.05 {(q<.05).sum()} | tasks {ok.sum()} | {time.time()-t0:.0f}s", flush=True)
