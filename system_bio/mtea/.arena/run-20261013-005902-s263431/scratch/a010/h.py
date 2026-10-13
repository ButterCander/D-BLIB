import pandas as pd, numpy as np, pickle, os
from scipy import stats
B = "D:/Cander/School/3-2/system_bio/mtea"
def P(*a): return os.path.join(B,*a)
auth = pickle.load(open("auth_moesm4.pkl","rb"))
def a_scores(sh):
    d = auth[sh].copy()
    d['tid'] = 'X'+d['Task ID'].astype(float).astype(int).astype(str)
    d = d.rename(columns={'TIDE score':'score','p-value':'p'})
    d['score']=d['score'].astype(float); d['p']=d['p'].astype(float)
    return d.set_index('tid')[['Description','score','p']]
def my(suf, tag="p10000"):
    return pd.read_csv(P("out",f"tide_{tag}_{suf}.tsv"), sep="\t").set_index('task_id')
def bh(p):
    p = np.asarray(p, float); n=len(p); o=np.argsort(p); r=np.empty(n); 
    q = p[o]*n/np.arange(1,n+1); q = np.minimum.accumulate(q[::-1])[::-1]; r[o]=np.minimum(q,1); return r
