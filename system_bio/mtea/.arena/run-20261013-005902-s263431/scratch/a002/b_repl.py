import numpy as np, pandas as pd
from scipy import stats
def load(p): return pd.read_csv(p,sep='\t').set_index('task_id')['score']
GS={"PI":"PI103","PD":"PD0325901","OXO":"5Z","PIOXO":"PI103+5Z","PIPD":"PI103+PD"}
T={k:load(f"out/tidex_tsir_{k}_24h.tsv") for k in GS}
G={k:load(f"out/tide_p10000_{v}.tsv") for k,v in GS.items()}
M=pd.DataFrame({c:{r:stats.spearmanr(*[pd.concat([G[r],T[c]],axis=1).dropna().iloc[:,i] for i in (0,1)])[0] for r in GS} for c in GS})
print(M.round(2))
d=np.diag(M.values); off=M.values[~np.eye(5,dtype=bool)]
print('matched',d.round(2),'mean',d.mean().round(3),'unmatched mean',off.mean().round(3),'max',off.max().round(2),'min',off.min().round(2))
# also with authors' scores vs Tsirvouli? 
