import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
M4='si/MOESM4.xlsx'
d={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    raw=pd.read_excel(M4,sheet_name=s,header=None)
    hdr=raw.iloc[0,2]
    df=raw.iloc[2:,:5].copy(); df.columns=['id','desc','score','p','sig']
    df=df.dropna(subset=['id'])
    df['score']=pd.to_numeric(df.score); df['p']=pd.to_numeric(df.p)
    d[s]=df
    p=df.p.values
    q=bh(np.clip(2*p,1e-12,1))
    print(s,'| header cell:',hdr,'| rows',len(df),'| p<0.025:',(df.p<0.025).sum(),'| sig flag:',df.sig.sum(),'| BH q<.05 (2p):',(q<.05).sum(),'| q<.10:',(q<.10).sum(),'| up',(df.score[df.p<0.025]>0).sum(),'down',(df.score[df.p<0.025]<0).sum())
# MOESM7
m7={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-MEKi','PI3Ki-TAKi']:
    x=pd.read_excel('si/MOESM7.xlsx',sheet_name=s)
    m7[s]=x
    print('M7',s,x.shape, list(x.columns)[:5])
# compare MOESM4 sheet scores vs MOESM7 tide scores
for s in ['TAKi','MEKi']:
    x=m7[s]; col=[c for c in x.columns if c.startswith('tide_') and not c.startswith('tide_e')][0]
    x=x.dropna(subset=['task_id']); x['id']=x.task_id.str.replace('X','').astype(float)
    for s4 in ['TAKi','MEKi']:
        a=d[s4]; j=a.merge(x,on='id'); 
        print('M7',s,'vs M4 sheet',s4,'n',len(j),'r',np.corrcoef(j.score.astype(float),pd.to_numeric(j[col],errors='coerce').fillna(0))[0,1])
# synergy
for s in ['PI3Ki-TAKi','PI3Ki-MEKi']:
    x=pd.read_excel('si/MOESM5.xlsx',sheet_name=s); x=x.dropna(subset=['task_id'])
    print('MOESM5',s,len(x), x.direction.value_counts().to_dict() if 'direction' in x else '')
print(pd.read_excel('si/MOESM5.xlsx',sheet_name='Summary',header=None).to_string())
