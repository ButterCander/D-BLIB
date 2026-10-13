import pandas as pd, numpy as np
from scipy import stats
B="D:/Cander/School/3-2/system_bio/mtea/"
def load(p): return pd.read_csv(B+p, sep="\t").set_index("task_id")["score"]
GS = {"PI": "PI103", "PD": "PD0325901", "OXO": "5Z", "PIOXO": "PI103+5Z", "PIPD": "PI103+PD"}
T = {k: load(f"out/tidex_tsir_{k}_24h.tsv") for k in GS}
G = {k: load(f"out/tide_p10000_{v}.tsv") for k, v in GS.items()}
TE = {k: load(f"out/tideex_tsir_{k}_24h.tsv") for k in GS}
GE = {k: load(f"out/tidee_{v}.tsv") for k, v in GS.items()}
def rho(a,b):
    j=pd.concat([a,b],axis=1).dropna(); return stats.spearmanr(j.iloc[:,0], j.iloc[:,1])[0]
for lab,X,Y in (("TIDE",G,T),("TIDE-ess",GE,TE)):
    M = pd.DataFrame({c:{r:rho(X[r],Y[c]) for r in GS} for c in GS})
    print(lab); print(M.round(2).to_string())
    d=np.diag(M.values); off=M.values[~np.eye(5,dtype=bool)]
    print("matched mean %.2f min %.2f max %.2f | unmatched mean %.2f max %.2f"%(d.mean(),d.min(),d.max(),off.mean(),off.max()))
