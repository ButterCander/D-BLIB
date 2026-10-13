import os, time, glob
import numpy as np, pandas as pd
from cobra.io import read_sbml_model
from mteapy.utils import mask_lfc_values
from mteapy.tide import compute_TIDE, compute_TIDEe
D = "mteapy/src/mteapy/data"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t"); model_ids = set(gm.geneID)
ts = pd.read_csv(f"{D}/task_structure_matrix.tsv", sep="\t", index_col=0)
ge = pd.read_csv(f"{D}/HumanGEM_essential_genes_matrix.tsv", sep="\t", index_col=0)
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
for f in sorted(glob.glob("out/statin_de_*.tsv")):
    name = os.path.basename(f)[len("statin_de_"):-4]
    d = pd.read_csv(f, sep="\t").rename(columns={"gene": "gene", "logFC": "lfc", "adj.P.Val": "padj"})
    d = d[d.gene.isin(model_ids)].drop_duplicates("gene")[["gene", "lfc", "padj"]]
    e = mask_lfc_values(d.copy(), "lfc", "padj", 0.05).set_index("gene")
    r = compute_TIDE(e, "lfc", ts, model, "absmax", 1000, 1, False, 1); r.to_csv(f"out/tidex_statin_{name}.tsv", sep="\t", index=False)
    re_ = compute_TIDEe(e, "lfc", ge, 1000, 1, False, 1); re_.to_csv(f"out/tideex_statin_{name}.tsv", sep="\t", index=False)
    print(f"{name}: model genes {len(e)}, nonzero {(e.lfc!=0).sum()}", flush=True)
