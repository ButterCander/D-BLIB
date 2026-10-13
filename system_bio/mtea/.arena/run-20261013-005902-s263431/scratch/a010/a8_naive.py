import json, numpy as np, pandas as pd
from scipy import stats
B="D:/Cander/School/3-2/system_bio/mtea/"
tasks = json.load(open(B+"out/task_genes.json"))
for c in ["PI103","PD0325901","5Z","PI103+5Z","PI103+PD"]:
    de = pd.read_csv(B+f"out/de_{c}.tsv", sep="\t", index_col=0)
    lf = de["lfc_shrunk"].where(de["padj"].fillna(1) < .05, 0.0)
    naive = pd.Series({t: lf.reindex(g).dropna().mean() for t, g in tasks.items() if lf.reindex(g).dropna().size >= 1})
    tide = pd.read_csv(B+f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")["score"]
    j = pd.concat([tide, naive], axis=1, keys=["tide","naive"]).dropna()
    print(c, len(j), "pearson %.2f spearman %.2f"%(stats.pearsonr(j.tide,j.naive)[0], stats.spearmanr(j.tide,j.naive)[0]))
