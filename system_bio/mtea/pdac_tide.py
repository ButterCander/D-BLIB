"""TIDE / TIDE-essential for (a) PDAC trametinib time course (GSE238198, limma-trend logFC) and (b) Tsirvouli 24 h official limma logFC.
Same pipeline as repro.py: mask logFC with padj >= 0.05 -> 0, model genes only, absmax GPR, NPERM label permutations. Run with the project venv."""
import os, sys, time
import numpy as np, pandas as pd
from cobra.io import read_sbml_model
from mteapy.utils import mask_lfc_values
from mteapy.tide import compute_TIDE, compute_TIDEe
NPERM = int(os.environ.get("NPERM", "2000"))
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
model_ids = set(gm.values)
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
ge = pd.read_csv(f"{D}/HumanGEM_essential_genes_matrix.tsv", sep="\t", index_col=0)
t0 = time.time(); model = read_sbml_model(f"{D}/HumanGEM.xml.gz"); print(f"model loaded {time.time()-t0:.0f}s", flush=True)

jobs = []   # (name, DataFrame[gene, lfc, padj])
for ln in ("CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"):
    for k in ("h3", "h12", "h24"):
        d = pd.read_csv(f"out/pdac_de_{ln}_{k}.tsv", sep="\t")
        d["gene"] = d.gene.map(gm); d = d.dropna(subset=["gene"]).drop_duplicates("gene")
        jobs.append((f"pdac_{ln}_{k}", d.rename(columns={"logFC": "lfc", "adj.P.Val": "padj"})[["gene", "lfc", "padj"]]))
off = pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv")
for c in ("PI", "PD", "OXO", "PIPD", "PIOXO"):
    d = off[off.Contrast == f"Diff_{c}_24h"].rename(columns={"Gene_name": "gene", "logFC": "lfc", "adj.P.Val": "padj"})
    d = d[d.gene.isin(model_ids)].drop_duplicates("gene")[["gene", "lfc", "padj"]]
    jobs.append((f"tsir_{c}_24h", d))
for name, d in jobs:
    e = mask_lfc_values(d.copy(), "lfc", "padj", 0.05).set_index("gene")
    t1 = time.time()
    r = compute_TIDE(e, "lfc", ts, model, "absmax", NPERM, 1, False, 1); r.to_csv(f"out/tidex_{name}.tsv", sep="\t", index=False)
    re_ = compute_TIDEe(e, "lfc", ge, NPERM, 1, False, 1); re_.to_csv(f"out/tideex_{name}.tsv", sep="\t", index=False)
    print(f"{name}: model genes {len(e)}, nonzero {(e.lfc!=0).sum()}, TIDE p<.025 {(r.pvalue<.025).sum()}, TIDE-e p<.025 {(re_.pvalue<.025).sum()} ({time.time()-t1:.0f}s)", flush=True)
