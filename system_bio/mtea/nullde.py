"""Negative control for the whole pipeline: shuffle condition labels within each batch (breaks all drug effects),
rerun DESeq2 (pydeseq2) -> task scores (TIDE, TIDE-essential) -> p-values. Expect ~5% raw hits and ~0 after BH."""
import os, sys, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, "mteapy", "src", "mteapy", "data"); OUT = os.path.join(HERE, "out")
SEED = int(os.environ.get("SEED", "1")); NPERM = int(os.environ.get("NPERM", "3000"))
CONDS = ["control", "PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]
rng = np.random.default_rng(SEED)
labels = []
for b in range(4):
    labels += list(rng.permutation(CONDS))          # each batch keeps one sample of every condition
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from mteapy.utils import mask_lfc_values
from mteapy.tide import compute_TIDE, compute_TIDEe
from cobra.io import read_sbml_model
cnt = pd.read_csv(os.path.join(HERE, "GSE285616_counts.csv.gz"), index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]
meta = pd.DataFrame({"batch": [f"b{i // 6 + 1}" for i in range(24)], "condition": pd.Categorical(labels, categories=CONDS)}, index=list(cnt.columns))
cnt = cnt.loc[cnt.sum(axis=1) >= 10].T
dds = DeseqDataSet(counts=cnt, metadata=meta, design="~batch + condition", refit_cooks=True, n_cpus=1, quiet=True); dds.deseq2()
gm = pd.read_csv(os.path.join(D, "HumanGEM_genes.tsv"), sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts = pd.read_csv(os.path.join(D, "task_structure_matrix.tsv"), sep="\t", index_col=0)
ge = pd.read_csv(os.path.join(D, "HumanGEM_essential_genes_matrix.tsv"), sep="\t", index_col=0)
model = read_sbml_model(os.path.join(D, "HumanGEM.xml.gz"))
for t in CONDS[1:]:
    st = DeseqStats(dds, contrast=["condition", t, "control"], n_cpus=1, quiet=True); st.summary()
    r = st.results_df.copy(); r.index = [gm.get(i, i) for i in r.index]; r.index.name = "gene"
    e = r.reset_index()[["gene", "log2FoldChange", "padj"]]; e["padj"] = e["padj"].fillna(1)
    e = mask_lfc_values(e, "log2FoldChange", "padj", 0.05).set_index("gene")
    ndeg = int((r.padj < .05).sum())
    a = compute_TIDEe(e, "log2FoldChange", ge, NPERM, 1, False, 1); b = compute_TIDE(e, "log2FoldChange", ts, model, "absmax", NPERM, 1, False, 1)
    a.to_csv(os.path.join(OUT, f"null{SEED}_tidee_{t}.tsv"), sep="\t", index=False); b.to_csv(os.path.join(OUT, f"null{SEED}_tide_{t}.tsv"), sep="\t", index=False)
    print(f"null seed {SEED} {t}: padj<.05 genes {ndeg} | TIDE-e p<.025 {int((a.pvalue<.025).sum())}/{len(a)} | TIDE p<.025 {int((b.pvalue<.025).sum())}/{len(b)}", flush=True)
