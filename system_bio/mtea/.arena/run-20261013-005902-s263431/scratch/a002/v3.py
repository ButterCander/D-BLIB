import pandas as pd, numpy as np, pickle
from scipy import stats
A=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a002/A.pkl','rb'))
M=pd.DataFrame({k:v.set_index('id').score for k,v in A.items()})
P=pd.DataFrame({k:v.set_index('id').p for k,v in A.items()})
print(M.shape)
print('Pearson'); print(M.corr().round(2))
print('Spearman'); print(M.corr(method='spearman').round(2))
# top tasks
print(M['PI3Ki'].sort_values().head(5))
print(A['PI3Ki'][A['PI3Ki'].id==76])
for k in A: print(k, A[k][A[k].id==76][['score','p']].values)
# fraction down among significant
for k in A:
    d=A[k]; s=d[d.p<.025]; print(k,len(s),'down frac',(s.score<0).mean().round(3))
# range of significant scores
for k in ('PI3Ki-MEKi','PI3Ki-TAKi'):
    d=A[k]; s=d[d.p<.025]; print(k,'abs score median sig',s.score.abs().median(),'max',s.score.abs().max())
