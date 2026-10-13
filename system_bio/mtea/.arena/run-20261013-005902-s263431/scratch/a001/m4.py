import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
SI='si/MOESM4.xlsx'
res={}
for sh in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel(SI,sheet_name=sh,header=None)
    hdr=d.iloc[0].tolist()
    d=d.iloc[2:,:4]; d.columns=['id','desc','score','p']
    d=d.dropna(subset=['id'])
    d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    res[sh]=d
    n=len(d); sig=(d.p<0.025).sum()
    q=bh(np.clip(2*d.p.values,1e-12,1))
    q1=bh(np.clip(d.p.values,1e-12,1))
    print(sh,'hdr',hdr[:3],'rows',n,'p<.025',sig,'neg',((d.p<0.025)&(d.score<0)).sum(),'pos',((d.p<0.025)&(d.score>0)).sum(),'BH(2p)<.05',(q<.05).sum(),'BH(2p)<.1',(q<.10).sum(),'BH(p)<.05',(q1<.05).sum(), 'unique ids',d.id.nunique(),'min p',d.p.min())
