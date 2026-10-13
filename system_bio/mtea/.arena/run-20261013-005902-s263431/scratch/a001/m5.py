import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
exec(open('.arena/run-20261013-005902-s263431/scratch/a001/m4b.py').read().split("for sh in")[0])
S=pd.read_excel('si/MOESM5.xlsx',sheet_name='Summary',header=None); print(S.to_string())
for combo,partner in [('PI3Ki-MEKi','MEKi'),('PI3Ki-TAKi','TAKi')]:
    d=pd.read_excel('si/MOESM5.xlsx',sheet_name=combo)
    d=d.dropna(subset=['task_id']); print(combo,'rows',len(d),d.direction.value_counts().to_dict())
    cA=m4(combo); pa=m4('PI3Ki'); 
    # note headers swapped: the sheet named partner (MEKi/TAKi) names are as sheet name
    pb=m4(partner)
    q=pd.Series(bh(np.clip(2*cA.p.values,1e-12,1)),index=cA.index)
    cols=list(d.columns)
    A,B,AB=cols[2],cols[3],cols[4]
    d['strong']=(d[A]==0)&(d[B]==0)
    d['margin']=np.where(d[AB]<0,d[AB]-np.minimum(d[A],d[B]),d[AB]-np.maximum(d[A],d[B]))
    print(' strong count',d.strong.sum(), 'weak',(~d.strong).sum())
    w=d[~d.strong]
    print(' weak margin median %.3f ; frac<0.1 %.2f ; frac<0.2 %.2f; min %.3f max %.3f'%(w.margin.abs().median(),(w.margin.abs()<0.1).mean(),(w.margin.abs()<0.2).mean(),w.margin.abs().min(),w.margin.abs().max()))
    # margin among tasks where both singles nonzero
    for _,r in d[d.strong].iterrows():
        t=r.task_id
        print('   strong',t,r.description[:40],'AB %.2f'%r[AB],'pPI3Ki %.3f'%pa.p[t],'p%s %.3f'%(partner,pb.p[t]),'pAB %.4f'%cA.p[t],'q %.3f'%q[t])
    print(' q<.05 among calls',(q[d.task_id]<.05).sum(),'q<.10',(q[d.task_id]<.10).sum())
