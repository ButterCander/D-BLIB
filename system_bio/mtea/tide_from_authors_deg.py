"""Score TIDE directly from the AUTHORS' Supplementary Data 1 DEG lists (|LFC|>=1, padj<.05; all other genes 0) and compare with their Supplementary Data 3."""
import numpy as np, pandas as pd
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
gids = set(pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").geneID)
meta = pd.read_csv(f"{D}/task_metadata.tsv", sep="\t").set_index("ID")
sheets = ["PI3Ki", "MEKi", "TAKi", "PI3Ki-TAKi", "PI3Ki-MEKi"]
for sh in sheets:
    deg = pd.read_excel("si/MOESM2.xlsx", sheet_name=sh)
    allg = sorted(gids)
    lfc = deg[deg.geneID.isin(gids)].drop_duplicates("geneID").set_index("geneID")["log2FoldChange"]
    e = pd.DataFrame({"lfc": 0.0}, index=allg); e.loc[lfc.index, "lfc"] = lfc.values
    r = compute_TIDE(e, "lfc", ts, model, "absmax", 20, 1, False, 1).set_index("task_id")["score"]
    a = pd.read_excel("si/MOESM4.xlsx", sheet_name=sh, header=None).iloc[2:, [0, 2]]; a.columns = ["id", "s"]
    a = pd.Series(pd.to_numeric(a.s).values, index=["X" + str(i) for i in a.id])
    j = pd.DataFrame({"auth": a, "deg": r}).dropna()
    big = (j.auth - j.deg).abs().nlargest(3)
    print(f"{sh:11} model genes among authors' DEGs: {len(lfc):4d} | r(authors' score, score from authors' DEG list) = {j.auth.corr(j.deg):.3f} | max |diff| {(j.auth-j.deg).abs().max():.2f} | exact matches (|diff|<1e-6): {int(((j.auth-j.deg).abs()<1e-6).sum())}/{len(j)}", flush=True)
    if sh == "MEKi": print("   largest diffs:", [(i, meta.loc[i, 'DESCRIPTION'][:28], round(j.auth[i], 2), round(j.deg[i], 2)) for i in big.index])
print("done")
