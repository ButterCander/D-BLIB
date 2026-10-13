"""Benedicto 2025 reproduction check: counts -> DESeq2 (pydeseq2, ~batch+condition) -> TIDE-essential / TIDE.
usage: python repro.py de | tidee | tide   (run from system_bio/mtea)"""
import gzip
import os
import sys
import time
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "mteapy", "src", "mteapy", "data")
OUT = os.path.join(HERE, "out")
os.makedirs(OUT, exist_ok=True)
CONDS = ["control", "PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]
TREATS = CONDS[1:]
NPERM, SEED = int(os.environ.get("NPERM", "1000")), 1
ALL = os.environ.get("ALLGENES") == "1"  # keep non-model genes in the permutation null
TAG = ("_all" if ALL else "") + (f"_p{NPERM}" if NPERM != 1000 else "")


def de():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    cnt = pd.read_csv(os.path.join(HERE, "GSE285616_counts.csv.gz"), index_col="Symbol")
    cnt = cnt[~cnt.index.duplicated()]
    cols = list(cnt.columns)
    meta = pd.DataFrame({"batch": [f"b{i // 6 + 1}" for i in range(24)], "condition": pd.Categorical([CONDS[i % 6] for i in range(24)], categories=CONDS)}, index=cols)
    cnt = cnt.loc[cnt.sum(axis=1) >= 10].T
    dds = DeseqDataSet(counts=cnt, metadata=meta, design="~batch + condition", refit_cooks=True, n_cpus=1, quiet=True)
    dds.deseq2()
    print("coefficients:", list(dds.obsm["design_matrix"].columns))
    gm = pd.read_csv(os.path.join(D, "HumanGEM_genes.tsv"), sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
    for t in TREATS:
        st = DeseqStats(dds, contrast=["condition", t, "control"], n_cpus=1, quiet=True)
        st.summary()
        r = st.results_df.copy()
        coef = f"condition[T.{t}]"
        try:
            st.lfc_shrink(coeff=coef)
            r["lfc_shrunk"] = st.results_df["log2FoldChange"]
        except Exception as e:
            print("  shrink failed:", t, repr(e)[:120])
            r["lfc_shrunk"] = r["log2FoldChange"]
        r.to_csv(os.path.join(OUT, f"de_{t}.tsv"), sep="\t")
        sig = ((r.padj < 0.05) & (r.log2FoldChange.abs() >= 1)).sum()
        print(f"{t}: {len(r)} genes, DEG(|lfc|>=1, padj<.05) = {sig}")


def load_de(t, lfc="lfc_shrunk"):
    r = pd.read_csv(os.path.join(OUT, f"de_{t}.tsv"), sep="\t", index_col=0)
    gm = pd.read_csv(os.path.join(D, "HumanGEM_genes.tsv"), sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
    if not ALL:
        r = r.loc[r.index.intersection(gm.index)]
    r.index = [gm.get(i, i) for i in r.index]
    r.index.name = "gene"
    r = r.reset_index()[["gene", lfc, "padj"]]
    r["padj"] = r["padj"].fillna(1)
    return r


def tide(kind):
    from mteapy.utils import mask_lfc_values, add_task_metadata
    from mteapy.tide import compute_TIDE, compute_TIDEe
    meta = pd.read_csv(os.path.join(D, "task_metadata.tsv"), sep="\t")
    if kind == "tidee":
        ge = pd.read_csv(os.path.join(D, "HumanGEM_essential_genes_matrix.tsv"), sep="\t", index_col=0)
    else:
        from cobra.io import read_sbml_model
        ts = pd.read_csv(os.path.join(D, "task_structure_matrix.tsv"), sep="\t", index_col=0)
        t0 = time.time()
        model = read_sbml_model(os.path.join(D, "HumanGEM.xml.gz"))
        print(f"model loaded {time.time() - t0:.0f}s: {len(model.reactions)} reactions")
    for t in TREATS:
        e = mask_lfc_values(load_de(t), "lfc_shrunk", "padj", 0.05).set_index("gene")
        t0 = time.time()
        if kind == "tidee":
            res = compute_TIDEe(e, "lfc_shrunk", ge, NPERM, 1, False, SEED)
        else:
            res = compute_TIDE(e, "lfc_shrunk", ts, model, "absmax", NPERM, 1, False, SEED)
        res.to_csv(os.path.join(OUT, f"{kind}{TAG}_{t}.tsv"), sep="\t", index=False)
        nsig = int((res["pvalue"] < 0.025).sum())
        print(f"{kind} {t}: {len(res)} tasks, p<0.025: {nsig}, {time.time() - t0:.0f}s", flush=True)


{"de": de, "tidee": lambda: tide("tidee"), "tide": lambda: tide("tide")}[sys.argv[1]]()
