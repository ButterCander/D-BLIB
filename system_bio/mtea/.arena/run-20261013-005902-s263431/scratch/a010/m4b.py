import sys, pickle, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
b='D:/Cander/School/3-2/system_bio/mtea/.arena/run-20261013-005902-s263431/scratch/a010/'
d=pickle.load(open(b+'m4.pkl','rb'))
for k,v in d.items():
    r=v[v.id==76]
    print(k, r[['score','p','sig']].values.tolist())
# BH
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); 
    r=p[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
print('BH on authors p (one-sided tail p as given):')
for k,v in d.items():
    q=bh(v.p.values); q2=bh(np.minimum(1,2*v.p.values))
    print(k,'p<.025',(v.p<0.025).sum(),'BH<.05',(q<0.05).sum(),'BH<.10',(q<0.10).sum(),'| two-sided 2p BH<.05',(q2<0.05).sum(),' p min',v.p.min(), 'n p==0',(v.p==0).sum())
