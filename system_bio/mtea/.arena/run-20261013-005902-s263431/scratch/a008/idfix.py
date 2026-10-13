import numpy as np, pandas as pd, sys
sys.path.insert(0,'mteapy/src')
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
COND = {"PI103": "PI3Ki", "PD0325901": "MEKi", "5Z": "TAKi", "PI103+5Z": "PI3Ki-TAKi", "PI103+PD": "PI3Ki-MEKi"}
def build(c, fixed):
    r0 = pd.read_csv(f"out/de_{c}.tsv", sep="\t", index_col=0); r0["padj"] = r0["padj"].fillna(1)
    ids = []
    for s in r0.index:
        if fixed and s in auth_map:
            g = [i for i in auth_map[s] if i in gem_ids]; ids.append(g[0] if g else None)
        else: ids.append(sym2gem.get(s))
    r0["id"] = ids; m = r0[r0.id.notna()].drop_duplicates("id")
    return pd.DataFrame({"lfc": np.where(m.padj < .05, m.lfc_shrunk, 0.0)}, index=m.id.values)
for c, sh in COND.items():
    a = pd.read_excel("si/MOESM4.xlsx", sheet_name=sh, header=None).iloc[2:, [0, 2, 3]]; a.columns = ["id", "s", "p"]
    a["tid"] = "X" + a.id.astype(str); a["s"] = pd.to_numeric(a.s)
    res=[]
    for fixed in (False, True):
        r = compute_TIDE(build(c, fixed), "lfc", ts, model, "absmax", 10, 1, False, 1)
        j = a.merge(r, left_on="tid", right_on="task_id")
        res.append(np.corrcoef(j.s, j.score)[0,1])
    print(f"{c:10} symbol-map r={res[0]:.3f}  ID-fix r={res[1]:.3f}", flush=True)
