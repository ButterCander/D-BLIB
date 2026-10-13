import pandas as pd, numpy as np, pickle
from scipy.stats import false_discovery_control as bh
A=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a002/A.pkl','rb'))
for combo,partner in (('PI3Ki-MEKi','MEKi'),('PI3Ki-TAKi','TAKi')):
    d=pd.read_excel('si/MOESM5.xlsx',sheet_name=combo)
    d=d.iloc[:,:6]; d.columns=['task_id','desc','A','B','AB','direction']
    d['tid']=d.task_id.str[1:].astype(int)
    print(combo,len(d)); 
    # classify strong: A==0 & B==0
    d['strong']=(d.A==0)&(d.B==0)
    print(' strong',d.strong.sum(),' weak',(~d.strong).sum())
    ab=A[combo].set_index('id'); a=A['PI3Ki'].set_index('id'); b=A[partner].set_index('id')
    d['pAB']=d.tid.map(ab.p); d['sAB']=d.tid.map(ab.score); d['pA']=d.tid.map(a.p); d['pB']=d.tid.map(b.p); d['sA']=d.tid.map(a.score); d['sB']=d.tid.map(b.score)
    # check consistency of MOESM5 with MOESM4
    print(' sAB match', np.allclose(d.AB, d.sAB, atol=1e-6))
    q=pd.Series(bh(np.clip(2*ab.p.values,1e-12,1)),index=ab.index)
    d['qAB']=d.tid.map(q)
    print(' BH q<.05 among calls',(d.qAB<.05).sum(),' q<.10',(d.qAB<.10).sum())
    st=d[d.strong]
    print(st[['task_id','desc','AB','pAB','pA','pB','qAB']].to_string())
    wk=d[~d.strong].copy()
    # margin over best single drug (same-direction): for down: min(A,B)-AB ; up: AB-max
    def marg(r):
        if r.direction.startswith('down'): return min(r.A,r.B)-r.AB
        return r.AB-max(r.A,r.B)
    wk['margin']=wk.apply(marg,axis=1)
    print(' weak n',len(wk),'median margin',wk.margin.median(),'frac<0.10',(wk.margin<0.10).mean(), 'frac<0.05',(wk.margin<0.05).mean())
    # how many weak calls have AB significant itself
    print(' weak with pAB<.025',(wk.pAB<.025).sum(),'of',len(wk), ' all calls pAB<.025',(d.pAB<.025).sum())
    # singles significant?
    print(' weak where best single drug is significant (p<.025):', ((wk.pA<.025)|(wk.pB<.025)).sum())
