"""Interaction (non-additive) contrasts on the same counts: AB - A - B (full-dose additive) and AB - .5A - .5B (dose-scaled).
usage: python inter.py de | tide | tidee   (NPERM, ALLGENES env as in repro.py)"""
import os, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "mteapy", "src", "mteapy", "data"); OUT = os.path.join(HERE, "out")
CONDS = ["control", "PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]
NPERM = int(os.environ.get("NPERM", "1000")); ALL = os.environ.get("ALLGENES") == "1"
# name -> {condition: weight} (control is the reference level, so its weight is implicit)
CONTR = {
    "PI3Ki_MEKi_add": {"PI103+PD": 1, "PI103": -1, "PD0325901": -1},
    "PI3Ki_MEKi_dose": {"PI103+PD": 1, "PI103": -.5, "PD0325901": -.5},
    "PI3Ki_TAKi_add": {"PI103+5Z": 1, "PI103": -1, "5Z": -1},
    "PI3Ki_TAKi_dose": {"PI103+5Z": 1, "PI103": -.5, "5Z": -.5},
}
HSA = {  # combination vs each single drug (highest-single-agent logic, no dose-shape assumption)
    "PI3Ki_MEKi_vsPI3Ki": {"PI103+PD": 1, "PI103": -1},
    "PI3Ki_MEKi_vsMEKi": {"PI103+PD": 1, "PD0325901": -1},
    "PI3Ki_TAKi_vsPI3Ki": {"PI103+5Z": 1, "PI103": -1},
    "PI3Ki_TAKi_vsTAKi": {"PI103+5Z": 1, "5Z": -1},
}
if os.environ.get("SET") == "hsa":
    CONTR = HSA

def de():
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats
    cnt = pd.read_csv(os.path.join(HERE, "GSE285616_counts.csv.gz"), index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]
    cols = list(cnt.columns)
    meta = pd.DataFrame({"batch": [f"b{i // 6 + 1}" for i in range(24)], "condition": pd.Categorical([CONDS[i % 6] for i in range(24)], categories=CONDS)}, index=cols)
    cnt = cnt.loc[cnt.sum(axis=1) >= 10].T
    dds = DeseqDataSet(counts=cnt, metadata=meta, design="~batch + condition", refit_cooks=True, n_cpus=1, quiet=True)
    dds.deseq2()
    names = list(dds.obsm["design_matrix"].columns); print("coefficients:", names, flush=True)
    for k, w in CONTR.items():
        v = np.zeros(len(names))
        for c, x in w.items(): v[names.index(f"condition[T.{c}]")] = x
        st = DeseqStats(dds, contrast=v, n_cpus=1, quiet=True); st.summary()
        r = st.results_df; r.to_csv(os.path.join(OUT, f"inter_{k}.tsv"), sep="\t")
        print(k, "genes", len(r), "padj<.05:", int((r.padj < .05).sum()), "|lfc|>=.5 & padj<.05:", int(((r.padj < .05) & (r.log2FoldChange.abs() >= .5)).sum()), flush=True)

def tide(kind):
    from mteapy.utils import mask_lfc_values
    from mteapy.tide import compute_TIDE, compute_TIDEe
    gm = pd.read_csv(os.path.join(D, "HumanGEM_genes.tsv"), sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
    if kind == "tidee":
        ge = pd.read_csv(os.path.join(D, "HumanGEM_essential_genes_matrix.tsv"), sep="\t", index_col=0)
    else:
        from cobra.io import read_sbml_model
        ts = pd.read_csv(os.path.join(D, "task_structure_matrix.tsv"), sep="\t", index_col=0)
        model = read_sbml_model(os.path.join(D, "HumanGEM.xml.gz"))
    for k in CONTR:
        r = pd.read_csv(os.path.join(OUT, f"inter_{k}.tsv"), sep="\t", index_col=0)
        if not ALL: r = r.loc[r.index.intersection(gm.index)]
        r.index = [gm.get(i, i) for i in r.index]; r.index.name = "gene"
        e = r.reset_index()[["gene", "log2FoldChange", "padj"]]; e["padj"] = e["padj"].fillna(1)
        e = mask_lfc_values(e, "log2FoldChange", "padj", 0.05).set_index("gene")
        t0 = time.time()
        res = compute_TIDEe(e, "log2FoldChange", ge, NPERM, 1, False, 1) if kind == "tidee" else compute_TIDE(e, "log2FoldChange", ts, model, "absmax", NPERM, 1, False, 1)
        tag = ("_all" if ALL else "") + (f"_p{NPERM}" if NPERM != 1000 else "")
        res.to_csv(os.path.join(OUT, f"inter_{kind}{tag}_{k}.tsv"), sep="\t", index=False)
        print(kind, k, "p<0.025:", int((res.pvalue < .025).sum()), f"{time.time()-t0:.0f}s", flush=True)

{"de": de, "tide": lambda: tide("tide"), "tidee": lambda: tide("tidee")}[sys.argv[1]]()
