import pandas as pd, numpy as np, pickle
from scipy.stats import pearsonr, spearmanr
B='D:/Cander/School/3-2/system_bio/mtea/'
res=pickle.load(open(B+'.arena/run-20261013-005902-s263431/scratch/a006/moesm4.pkl','rb'))
mp={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
A={k:res[k].set_index('id') for k in res}
M={k:pd.read_csv(B+'out/tide_%s.tsv'%v,sep='\t').assign(id=lambda d:d.task_id.str[1:].astype(int)).set_index('id') for k,v in mp.items()}
sig=[set(A[k].index[A[k].p<0.025]) for k in A]
shared=sorted(set.intersection(*sig))
print('tasks significant (p<.025) in all five conditions (authors):',len(shared),shared)
for t in shared: print('  X%d'%t, A['PI3Ki'].desc[t][:60])
xs=[];ys=[]
for t in shared:
    for k in mp:
        xs.append(float(A[k].score[t])); ys.append(float(M[k].score[t]))
print('r over these %d values = %.3f'%(len(xs),pearsonr(xs,ys)[0]))
# the five conditions r over all 187
for k in mp:
    j=A[k].index
    print(k,'r all tasks %.3f'%pearsonr(A[k].score.astype(float),M[k].score[j])[0],' spearman %.3f'%spearmanr(A[k].score.astype(float),M[k].score[j])[0])
print('--- leave-one-task-out r (6 tasks x 5 cond)')
import itertools
for drop in shared:
    ts=[t for t in shared if t!=drop]
    xs=[float(A[k].score[t]) for t in ts for k in mp]; ys=[float(M[k].score[t]) for t in ts for k in mp]
    print('drop X%d -> r=%.4f'%(drop,pearsonr(xs,ys)[0]))
# per-task values
for t in shared:
    print('X%d'%t,' '.join('%s %.2f/%.2f'%(k,float(A[k].score[t]),float(M[k].score[t])) for k in mp))
