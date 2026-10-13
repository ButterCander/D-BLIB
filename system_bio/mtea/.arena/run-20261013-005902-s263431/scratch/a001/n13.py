import pandas as pd, numpy as np
from scipy.stats import spearmanr, pearsonr
A={'PI':'PI103','PD':'PD0325901','OXO':'5Z','PIOXO':'PI103+5Z','PIPD':'PI103+PD'}
mine={k:pd.read_csv(f'out/tide_p10000_{v}.tsv',sep='\t').set_index('task_id').score for k,v in A.items()}
tsir={k:pd.read_csv(f'out/tidex_tsir_{k}_24h.tsv',sep='\t').set_index('task_id').score for k in A}
M=pd.DataFrame(np.zeros((5,5)),index=list(A),columns=list(A))
for i in A:
    for j in A:
        a,b=mine[i].align(tsir[j],join='inner'); M.loc[i,j]=spearmanr(a,b)[0]
print(M.round(2))
d=np.diag(M.values); off=M.values[~np.eye(5,dtype=bool)]
print('matched',d.round(2),'range %.2f-%.2f mean %.2f'%(d.min(),d.max(),d.mean()))
print('unmatched mean %.2f max %.2f min %.2f'%(off.mean(),off.max(),off.min()))
# TIDE-essential
Me=pd.DataFrame(np.zeros((5,5)),index=list(A),columns=list(A))
mine_e={k:pd.read_csv(f'out/tidee_p10000_{v}.tsv',sep='\t').set_index('task_id').score for k,v in A.items()}
tsir_e={k:pd.read_csv(f'out/tideex_tsir_{k}_24h.tsv',sep='\t').set_index('task_id').score for k in A}
for i in A:
    for j in A:
        a,b=mine_e[i].align(tsir_e[j],join='inner'); a=a.dropna(); b=b.reindex(a.index).dropna(); a=a.reindex(b.index); Me.loc[i,j]=spearmanr(a,b)[0]
print(Me.round(2)); d=np.diag(Me.values); off=Me.values[~np.eye(5,dtype=bool)]
print('E matched',d.round(2),'mean %.2f min %.2f; unmatched mean %.2f max %.2f'%(d.mean(),d.min(),off.mean(),off.max()))
