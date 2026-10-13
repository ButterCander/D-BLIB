import pandas as pd, numpy as np
from scipy import stats
B="D:/Cander/School/3-2/system_bio/mtea/"
f = pd.read_csv(B+"out/flobak_synergy.tsv", sep="\t")
syn = f[f.call=="syn"]; non = f[f.call=="non"]
print("n syn",len(syn),"n non",len(non))
for col in ["hsa_margin_full","bliss_excess","bliss_excess_dir"]:
    print(col, "syn range %.2f..%.2f"%(syn[col].min(),syn[col].max()), "non range %.2f..%.2f"%(non[col].min(),non[col].max()), "mean/sd non %.3f / %.3f"%(non[col].mean(), non[col].std(ddof=1)))
    from sklearn.metrics import roc_auc_score if False else None
