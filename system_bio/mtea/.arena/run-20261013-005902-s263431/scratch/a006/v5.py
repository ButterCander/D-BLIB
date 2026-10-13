import numpy as np, pandas as pd
from scipy import stats
model = pd.read_csv("external/Model_24Q4.csv").set_index("ModelID")
n2i = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
E = pd.read_csv("external/expr_pancreas.tsv", sep="\t", index_col=0)
TR = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
def base(k): return E.loc[n2i[k], TR].mean()
sets={
 'doc composite (n=7)':{"CAPAN2":1.06,"CFPAC1":(1.31+0.05)/2,"PANC1":(1.12+0.33)/2,"MIAPACA2":(0.33-0.02)/2,"KP4":0.18,"ASPC1":-0.10,"SW1990":-0.17},
 'trametinib only, mean over its datasets (n=5)':{"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":(0.33-0.02)/2,"KP4":0.18},
 'trametinib 24h GSE238198 only (n=4)':{"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":0.33},
 'KRAS inhibitor 24h GSE279749 only (n=4)':{"CFPAC1":0.05,"PANC1":0.33,"ASPC1":-0.10,"SW1990":-0.17},
}
for nm,s in sets.items():
    ks=list(s); x=[base(k) for k in ks]; y=[s[k] for k in ks]
    r,p=stats.spearmanr(x,y); print(f'{nm}: rho={r:+.2f} p={p:.3f}  baseline',{k:round(base(k),2) for k in ks})
# trametinib + KRASi both measured in 2 lines
print('PANC1 trametinib 1.12 vs MRTX1133 0.33; CFPAC1 trametinib 1.31 vs BI-2865 0.05')
