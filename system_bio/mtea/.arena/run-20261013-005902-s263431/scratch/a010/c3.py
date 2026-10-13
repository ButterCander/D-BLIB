import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
out={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=s,header=None).iloc[2:,:5]
    d.columns=['id','desc','score','p','sig']; d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    d['sig']=d.sig.astype(str)
    out[s]=d.set_index('id')
    print(s,'n',len(d),'p<0.025',(d.p<0.025).sum(),'sig=True',(d.sig=='True').sum(),'p<=0.025',(d.p<=0.025).sum(), 'neg among sig', ((d.p<0.025)&(d.score<0)).sum())
for t in [138,140,141]:
    for s,d in out.items():
        r=d.loc[t]; print(t,r.desc,s,round(r.score,3),r.p,r.sig)
