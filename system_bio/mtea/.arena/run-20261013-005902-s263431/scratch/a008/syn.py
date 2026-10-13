import pandas as pd, numpy as np, pickle
from scipy.stats import false_discovery_control as bh
m4=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a008/m4.pkl','rb'))
x=pd.ExcelFile('si/MOESM5.xlsx')
for sh,(a,b,ab) in {'PI3Ki-TAKi':('PI3Ki','TAKi','PI3Ki-TAKi'),'PI3Ki-MEKi':('PI3Ki','MEKi','PI3Ki-MEKi')}.items():
    d=pd.read_excel(x,sh).iloc[:,:6]; d.columns=['task','desc',a,b,ab,'dir']
    d['id']=d.task.str[1:].astype(int)
    pab=m4[ab].set_index('id'); 
    q=pd.Series(bh(np.clip(2*pab.p.values,1e-12,1)),index=pab.index)
    d['p_ab']=d.id.map(pab.p); d['q_ab']=d.id.map(q)
    d['p_a']=d.id.map(m4[a].set_index('id').p); d['p_b']=d.id.map(m4[b].set_index('id').p)
    d['strong']=(d[a]==0)&(d[b]==0)
    d['margin']=d[ab].abs()-d[[a,b]].abs().max(axis=1)
    w=d[~d.strong]
    print(sh,'n',len(d),'strong',d.strong.sum(),'weak',len(w),'| weak margin median',round(w.margin.median(),3),'frac<0.10',round((w.margin<0.10).mean(),2),'frac<0.05',round((w.margin<0.05).mean(),2),'| q_ab<.05',int((d.q_ab<.05).sum()),' q<.10',int((d.q_ab<.10).sum()))
    print(d[d.strong][['task','desc','p_a','p_b','p_ab','q_ab']].round(3).to_string())
    print('  weak: ab p<.025 count', int((w.p_ab<.025).sum()),'of',len(w), '; weak where single drug p<.025 (the non-zero one):', int(((w.p_a<.025)|(w.p_b<.025)).sum()))
