import numpy as np, pandas as pd
from scipy import stats
P = r"D:/Cander/School/3-2/system_bio/mtea/"
model = pd.read_csv(P + "external/Model_24Q4.csv").set_index("ModelID")
n2i = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
E = pd.read_csv(P + "external/expr_pancreas.tsv", sep="\t", index_col=0)
TR = ["HMGCS1", "HMGCR", "MVK", "MVD", "IDI1", "FDPS"]
lines = ["CAPAN2", "CFPAC1", "PANC1", "MIAPACA2", "KP4", "ASPC1", "SW1990"]
base = {k: E.loc[n2i[k], TR].mean() for k in lines}
print({k: round(v, 2) for k, v in base.items()})
IND = {"CAPAN2": [1.057], "CFPAC1": [1.313, 0.05], "PANC1": [1.125, 0.33], "MIAPACA2": [0.327, -0.02], "KP4": [0.18], "ASPC1": [-0.10], "SW1990": [-0.17]}
def rho(ind, ks):
    return stats.spearmanr([base[k] for k in ks], [ind[k] for k in ks])
ind = {k: np.mean(v) for k, v in IND.items()}
print('doc composite n=7', rho(ind, lines))
# trametinib only (24h GSE238198 for 4 lines; KP4 48h)
tram24 = {"CAPAN2": 1.057, "CFPAC1": 1.313, "PANC1": 1.125, "MIAPACA2": 0.327}
print('GSE238198 only (n=4)', rho(tram24, list(tram24)))
tram_all = dict(tram24); tram_all["KP4"] = 0.18
print('trametinib datasets, KP4 added (n=5)', rho(tram_all, list(tram_all)))
# KRAS inhibitors only (GSE279749)
kras = {"CFPAC1": 0.05, "PANC1": 0.33, "ASPC1": -0.10, "SW1990": -0.17}
print('GSE279749 only (n=4)', rho(kras, list(kras)))
# use the MAPK-pathway max or the first-listed dataset per line (a different plausible composite)
alt = {"CAPAN2": 1.057, "CFPAC1": 1.313, "PANC1": 1.125, "MIAPACA2": 0.327, "KP4": 0.18, "ASPC1": -0.10, "SW1990": -0.17}
print('trametinib 24h value where available else other (n=7)', rho(alt, lines))
# Exact permutation p for rho=0.964, n=7
import itertools
x = np.arange(7); cnt = 0; tot = 0; obs = stats.spearmanr(x, [0,1,2,3,4,6,5])[0]
for p in itertools.permutations(range(7)):
    r = stats.spearmanr(x, p)[0]; tot += 1
    if abs(r) >= 0.9643: cnt += 1
print('exact two-sided p for |rho|>=0.964, n=7:', cnt / tot, cnt, tot)
# Stability: baseline in GSE238198 NT samples
TR2 = TR
cnt_ = {}
for l in ("CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"):
    d = pd.read_csv(P + f"pdac/GSE238198_Pdac.{l}.genes.normalizedCounts.txt.gz", sep="\t", index_col=0)
    cnt_[l] = np.log2(d[[c for c in d.columns if "_NT_" in c]] + 1).mean(axis=1).reindex(TR2).mean()
print({k: round(v, 2) for k, v in cnt_.items()})
