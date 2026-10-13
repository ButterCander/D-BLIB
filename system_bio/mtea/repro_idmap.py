"""Fix: map DE results to Human-GEM by the authors' Ensembl IDs (Supplementary Data 1) where available, instead of by gene symbol."""
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
hg = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol")
sym2gem = hg.set_index("geneSymbol")["geneID"]; gem_ids = set(hg.geneID)
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
x = pd.ExcelFile("si/MOESM2.xlsx")
auth_map = {}
for sh in x.sheet_names:
    d = x.parse(sh)
    for s, g in zip(d.geneSymbol, d.geneID): auth_map.setdefault(s, set()).add(g)
conflict = {s: g for s, g in auth_map.items() if s in sym2gem.index and sym2gem[s] not in g}
print(f"symbols in authors' DEG lists: {len(auth_map)}; whose ID differs from Human-GEM's ID for that symbol: {len(conflict)}", flush=True)
print("  examples:", list(conflict.items())[:6], flush=True)
COND = {"PI103": "PI3Ki", "PD0325901": "MEKi", "5Z": "TAKi", "PI103+5Z": "PI3Ki-TAKi", "PI103+PD": "PI3Ki-MEKi"}
def build(c, fixed, N=None):
    r0 = pd.read_csv(f"out/de_{c}.tsv", sep="\t", index_col=0); r0["padj"] = r0["padj"].fillna(1)
    ids = []
    for s in r0.index:
        if fixed and s in auth_map:
            g = [i for i in auth_map[s] if i in gem_ids]            # authors' ID(s) that exist in the model; none -> gene not a model gene
            ids.append(g[0] if g else None)
        else:
            ids.append(sym2gem.get(s))
    r0["id"] = ids; m = r0[r0.id.notna()].drop_duplicates("id")
    e = pd.DataFrame({"lfc": np.where(m.padj < .05, m.lfc_shrunk, 0.0)}, index=m.id.values)
    if N and N > len(e): e = pd.concat([e, pd.DataFrame({"lfc": np.zeros(N - len(e))}, index=[f"d{i}" for i in range(N - len(e))])])
    return e
print("\ncondition | r(scores) symbol-map -> ID-fix | p<.025 / BH at N=5,000 zero-padded (authors' own in brackets)", flush=True)
for c, sh in COND.items():
    a = pd.read_excel("si/MOESM4.xlsx", sheet_name=sh, header=None).iloc[2:, [0, 2, 3]]; a.columns = ["id", "s", "p"]
    a["tid"] = "X" + a.id.astype(str); a["s"] = pd.to_numeric(a.s); a["p"] = pd.to_numeric(a.p)
    out = []
    for fixed in (False, True):
        e = build(c, fixed, 5000)
        r = compute_TIDE(e, "lfc", ts, model, "absmax", 3000, 1, False, 1)
        j = a.merge(r, left_on="tid", right_on="task_id")
        q = bh(np.clip(2 * j.pvalue.values, 1e-12, 1))
        out.append((np.corrcoef(j.s, j.score)[0, 1], int((j.pvalue < .025).sum()), int((q < .05).sum())))
    qa = bh(np.clip(2 * a.p.values, 1e-12, 1))
    print(f"{c:10} r {out[0][0]:.2f} -> {out[1][0]:.2f} | symbol-map {out[0][1]}/{out[0][2]} -> ID-fix {out[1][1]}/{out[1][2]}  [authors {(a.p<.025).sum()}/{(qa<.05).sum()}]", flush=True)
print("done")
