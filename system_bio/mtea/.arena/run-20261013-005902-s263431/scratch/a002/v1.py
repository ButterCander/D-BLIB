import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
sheets=['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']
A={}
for s in sheets:
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=s,header=None,skiprows=2).iloc[:,:4]
    d.columns=['id','desc','score','p']
    d=d[pd.to_numeric(d.id,errors='coerce').notna()].copy()
    d['id']=d.id.astype(float).astype(int); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    A[s]=d
    n=len(d); sig=(d.p<0.025).sum(); 
    q2=bh(np.clip(2*d.p.values,1e-12,1)); q1=bh(np.clip(d.p.values,1e-12,1))
    print(s,'n',n,'p<.025',sig,'up',((d.p<.025)&(d.score>0)).sum(),'down',((d.p<.025)&(d.score<0)).sum(),'BH2p<.05',(q2<.05).sum(),'BH2p<.10',(q2<.10).sum(),'BH(raw)<.05',(q1<.05).sum(),'p==0',(d.p==0).sum(), 'minp',d.p.min())
import pickle; pickle.dump(A,open('.arena/run-20261013-005902-s263431/scratch/a002/A.pkl','wb'))
