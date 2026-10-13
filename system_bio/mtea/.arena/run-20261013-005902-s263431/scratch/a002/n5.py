import pandas as pd, numpy as np, pickle
from scipy.stats import false_discovery_control as bh
A=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a002/A.pkl','rb'))
M={k:v.set_index('id').score for k,v in A.items()}
mp={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
for tag in ['tide_p10000','tide_all_p10000']:
    print(tag)
    for sh,f in mp.items():
        d=pd.read_csv(f'out/{tag}_{f}.tsv',sep='\t'); d['id']=d.task_id.str[1:].astype(int); d=d.set_index('id')
        row=[]
        for sh2 in M:
            j=pd.concat([M[sh2],d.score],axis=1,keys=['a','m']).dropna()
            row.append(f"{sh2}:{np.corrcoef(j.a,j.m)[0,1]:.2f}")
        q=bh(np.clip(2*d.pvalue.values,1e-12,1))
        print(f'  mine {f:10}', ' '.join(row), '| p<.025',(d.pvalue<.025).sum(),'BH',(q<.05).sum())
