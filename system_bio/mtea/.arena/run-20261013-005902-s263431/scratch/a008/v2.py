import pandas as pd, numpy as np, pickle
from scipy.stats import false_discovery_control as bh
M=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a008/M4.pkl','rb'))
x=pd.ExcelFile('si/MOESM5.xlsx')
print(pd.read_excel(x,'Summary',header=None).to_string())
for sh,(a,b) in {'PI3Ki-TAKi':('PI3Ki','TAKi'),'PI3Ki-MEKi':('PI3Ki','MEKi')}.items():
    d=pd.read_excel(x,sh); d=d.iloc[:,:6]; d.columns=['task','desc',a,b,'ab','dir']
    d=d.dropna(subset=['task']); d=d[d.task.astype(str).str.startswith('X')]
    d['id']=d.task.str[1:].astype(int)
    print(sh,'rows',len(d), 'dir counts', d.dir.value_counts().to_dict())
    pab=M[sh].set_index('id'); q=pd.Series(bh(np.clip(2*pab.p.values,1e-12,1)),index=pab.index)
    d['p_ab']=d.id.map(pab.p); d['q_ab']=d.id.map(q)
    d['pa']=d.id.map(M[a].set_index('id').p); d['pb']=d.id.map(M[b].set_index('id').p)
    strong=(d[a]==0)&(d[b]==0)
    print(' strong',strong.sum(),'weak',(~strong).sum(),' q<.05',(d.q_ab<.05).sum(),' q<.10',(d.q_ab<.10).sum(),' p_ab<.025',(d.p_ab<.025).sum())
    w=d[~strong].copy(); w['margin']=w.ab.abs()-w[[a,b]].abs().max(axis=1)
    print(' weak margin median %.3f  frac<0.10 %.2f  frac<0.05 %.2f' % (w.margin.median(), (w.margin<.10).mean(), (w.margin<.05).mean()))
    s=d[strong]
    print(s[['task','desc','pa','pb','p_ab','q_ab']].round(3).to_string())
    # of the weak ones: how many have the combination p<0.025?
    # strong: single drug p between .025 and .10?
    print(' strong single p min over both:', s[['pa','pb']].min(axis=1).round(3).tolist())
