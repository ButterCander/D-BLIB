"""Masking-rule variants for the PI103 reproduction (model-gene null, 3000 permutations): what narrows the null to the authors' width?"""
import time, numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
r0 = pd.read_csv("out/de_PI103.tsv", sep="\t", index_col=0)
r0 = r0.loc[r0.index.intersection(gm.index)]; r0.index = [gm[i] for i in r0.index]; r0.index.name = "gene"; r0["padj"] = r0["padj"].fillna(1)
au = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki", header=None).iloc[2:, [0, 2, 3]]
au.columns = ["tid", "sa", "pa"]; au["tid"] = au.tid.astype(str); au["pa"] = pd.to_numeric(au.pa); au["sa"] = pd.to_numeric(au.sa)
sd = lambda x, p: np.abs(x) / stats.norm.isf(np.clip(p, 1e-4, 0.5))
V = {
 "padj<.05 (baseline)":            lambda d: (d.padj < .05),
 "padj<.05 & |shrunk|>=1":         lambda d: (d.padj < .05) & (d.lfc_shrunk.abs() >= 1),
 "padj<.05 & |shrunk|>=0.5":       lambda d: (d.padj < .05) & (d.lfc_shrunk.abs() >= .5),
 "padj<.01":                       lambda d: (d.padj < .01),
 "padj<.05 & |unshrunk|>=1":       lambda d: (d.padj < .05) & (d.log2FoldChange.abs() >= 1),
}
print(f"authors: p<.025 {(au.pa<.025).sum()}, implied null sd {np.nanmedian(sd(au.sa, au.pa)):.3f}, sd of scores {au.sa.std():.3f}")
for name, rule in V.items():
    e = r0.reset_index()[["gene", "lfc_shrunk", "padj"]].copy()
    keep = rule(r0).reindex(r0.index).values
    e["lfc_shrunk"] = np.where(keep, e.lfc_shrunk, 0.0); e = e.set_index("gene")
    t0 = time.time(); res = compute_TIDE(e, "lfc_shrunk", ts, model, "absmax", 3000, 1, False, 1)
    res["tid"] = res.task_id.astype(str).str.lstrip("X"); j = au.merge(res, on="tid")
    q = bh(np.clip(2 * j.pvalue.values, 1e-12, 1)); slope = np.polyfit(j.score, j.sa, 1)[0]
    print(f"{name:26} nonzero {(e.lfc_shrunk!=0).sum():5d} | r {np.corrcoef(j.sa, j.score)[0,1]:.2f} slope(auth~mine) {slope:.2f} | p<.025 {(j.pvalue<.025).sum():3d} BH {(q<.05).sum():3d} | null sd {np.nanmedian(sd(j.score, j.pvalue)):.3f} | score sd {j.score.std():.3f} | {time.time()-t0:.0f}s", flush=True)
print("done")
