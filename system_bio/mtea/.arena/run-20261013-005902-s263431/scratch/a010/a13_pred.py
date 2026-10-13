import numpy as np, pandas as pd
from scipy import stats
B="D:/Cander/School/3-2/system_bio/mtea/"
model = pd.read_csv(B+"external/Model_24Q4.csv").set_index("ModelID")
name2id = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
E = pd.read_csv(B+"external/expr_pancreas.tsv", sep="\t", index_col=0)
IND = {"CAPAN2": [1.06], "CFPAC1": [1.31, 0.05], "PANC1": [1.12, 0.33], "MIAPACA2": [0.33, -0.02], "KP4": [0.18], "ASPC1": [-0.10], "SW1990": [-0.17]}
ind = pd.Series({k: np.mean(v) for k,v in IND.items()})
TR = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
ids = {k:name2id[k] for k in IND}
bt = pd.Series({k: E.loc[ids[k], TR].mean() for k in IND})
print(pd.DataFrame({"ind":ind,"baseline_trunk_log2TPM1":bt}).sort_values("ind",ascending=False).round(2))
r,p = stats.spearmanr(ind, bt); print("rho n=7: %.3f p=%.4f"%(r,p))
# exact permutation p for n=7
import itertools
cnt=0;tot=0
ri = stats.rankdata(ind.values); rb=stats.rankdata(bt.values)
for perm in itertools.permutations(range(7)):
    rr = np.corrcoef(ri, rb[list(perm)])[0,1]; tot+=1; cnt += abs(rr)>=abs(r)-1e-12
print("exact two-sided perm p = %.5f"%(cnt/tot))
# within GSE238198 only (n=4) – baseline from NT of the same study and from DepMap
ind238 = {"CAPAN2":1.06,"CFPAC1":1.31,"MIAPACA2":0.33,"PANC1":1.12}
x=[bt[k] for k in ind238]; y=list(ind238.values())
print("DepMap baseline vs 24h trametinib induction within GSE238198 (n=4): rho %.2f"%stats.spearmanr(x,y)[0])
# study composition of inducibility
print("\nWhich dataset contributes each line's inducibility:")
src = {"CAPAN2":"GSE238198 tram 24h only","CFPAC1":"GSE238198 tram 24h (+1.31) + GSE279749 BI-2865 (+0.05)","PANC1":"GSE238198 tram 24h (+1.12) + GSE279749 MRTX1133 (+0.33)","MIAPACA2":"GSE238198 tram 24h (+0.33) + GSE206267 tram 48h (-0.02)","KP4":"GSE206267 tram 48h only","ASPC1":"GSE279749 MRTX1133 only","SW1990":"GSE279749 MRTX1133 only"}
for k in ind.sort_values(ascending=False).index: print(f"  {k:9s} {ind[k]:+.2f}  {src[k]}")
# leave-out: using only trametinib datasets (GSE238198 24h + GSE206267 48h)
tram = {"CAPAN2":[1.06],"CFPAC1":[1.31],"PANC1":[1.12],"MIAPACA2":[0.33,-0.02],"KP4":[0.18]}
it = pd.Series({k:np.mean(v) for k,v in tram.items()})
print("\ntrametinib-only inducibility (5 lines):", it.round(2).to_dict())
print("rho baseline vs trametinib-only (n=5): %.2f"%stats.spearmanr(it, bt[it.index])[0])
# KRAS inhibitor-only lines (AsPC-1, SW1990, PANC-1, CFPAC-1) from GSE279749
