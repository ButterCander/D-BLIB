"""Why is my null wider than the authors'? PI103 only, 3000 permutations, model-gene null; variants of the LFC input."""
import os, time, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.utils import mask_lfc_values
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
r0 = pd.read_csv("out/de_PI103.tsv", sep="\t", index_col=0)
r0 = r0.loc[r0.index.intersection(gm.index)]; r0.index = [gm[i] for i in r0.index]; r0.index.name = "gene"
r0["padj"] = r0["padj"].fillna(1)
au = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki", header=None).iloc[2:, [0, 2, 3]]
au.columns = ["tid", "sa", "pa"]; au["tid"] = au.tid.astype(str); au["pa"] = pd.to_numeric(au.pa); au["sa"] = pd.to_numeric(au.sa)
sd = lambda x, p: np.abs(x) / stats.norm.isf(np.clip(p, 1e-4, 0.5))
print(f"genes {len(r0)}; |shrunk|>0 after mask {(mask_lfc_values(r0.reset_index()[['gene','lfc_shrunk','padj']], 'lfc_shrunk', 'padj', 0.05)['lfc_shrunk'] != 0).sum()}")
for name, col, thr in (("shrunk, mask padj<.05 (baseline)", "lfc_shrunk", 0.05), ("unshrunk, mask padj<.05", "log2FoldChange", 0.05), ("shrunk, no mask", "lfc_shrunk", 1.01)):
    e = r0.reset_index()[["gene", col, "padj"]]
    e = mask_lfc_values(e, col, "padj", thr).set_index("gene")
    t0 = time.time(); res = compute_TIDE(e, col, ts, model, "absmax", 3000, 1, False, 1)
    res["tid"] = res.task_id.astype(str).str.lstrip("X"); j = au.merge(res, on="tid")
    q = bh(np.clip(2 * j.pvalue, 1e-12, 1))
    print(f"{name:34} nonzero {(e[col]!=0).sum():5d} | r_vs_authors {np.corrcoef(j.sa, j.score)[0,1]:.2f} | p<.025 {(j.pvalue<.025).sum():3d} BH {(q<.05).sum():3d} | implied null sd {np.nanmedian(sd(j.score, j.pvalue)):.3f} (authors {np.nanmedian(sd(j.sa, j.pa)):.3f}) | {time.time()-t0:.0f}s", flush=True)
