import numpy as np, pandas as pd, sys
from scipy import stats
sys.stdout.reconfigure(encoding="utf-8")
model = pd.read_csv("external/Model_24Q4.csv").set_index("ModelID")
n2i = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
E = pd.read_csv("external/expr_pancreas.tsv", sep="\t", index_col=0)
TR = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
bt = {k: E.loc[n2i[k], TR].mean() for k in ["CAPAN2","CFPAC1","PANC1","MIAPACA2","KP4","ASPC1","SW1990"]}
print({k: round(v,2) for k,v in bt.items()})
comp = {"CAPAN2":1.06,"CFPAC1":0.68,"PANC1":0.725,"MIAPACA2":0.155,"KP4":0.18,"ASPC1":-0.10,"SW1990":-0.17}
def rho(d, label):
    ks=list(d); r,p=stats.spearmanr([bt[k] for k in ks],[d[k] for k in ks]); print(f"{label}: n={len(ks)} rho={r:+.2f} p={p:.4f}")
rho(comp,"composite (as in doc)")
tram24 = {"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":0.33}
rho(tram24,"trametinib 24h only (GSE238198)")
tram5 = dict(tram24); tram5["KP4"]=0.18
rho(tram5,"trametinib incl KP4 48h")
mrtx = {"PANC1":0.33,"ASPC1":-0.10,"SW1990":-0.17,"CFPAC1":0.05}
rho(mrtx,"KRAS-inhibitor 24h only (GSE279749)")
# 7 lines but with trametinib where available & KRAS-i otherwise (a different composite)
alt = {"CAPAN2":1.06,"CFPAC1":1.31,"PANC1":1.12,"MIAPACA2":0.33,"KP4":0.18,"ASPC1":-0.10,"SW1990":-0.17}
rho(alt,"tram where available, MRTX otherwise")
# only the two lines w/out trametinib removed => n=5 composite
rho({k:comp[k] for k in ["CAPAN2","CFPAC1","PANC1","MIAPACA2","KP4"]},"composite without AsPC-1/SW1990")
# exact p for n=7 rho 0.96
from itertools import permutations
x=np.arange(7); y=np.array(sorted(comp.values()))
