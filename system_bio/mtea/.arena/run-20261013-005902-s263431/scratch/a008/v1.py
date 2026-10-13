import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
x=pd.ExcelFile('si/MOESM4.xlsx'); M={}
for s in x.sheet_names:
    d=pd.read_excel(x,s,header=None)
    hdr=d.iloc[0,2]
    t=d.iloc[2:,:5].copy(); t.columns=['id','desc','score','p','sig']
    t['id']=pd.to_numeric(t.id,errors='coerce'); t=t.dropna(subset=['id']); t['id']=t.id.astype(int)
    t['score']=pd.to_numeric(t.score); t['p']=pd.to_numeric(t.p)
    M[s]=t
    n=(t.p<0.025).sum(); ns=int((t.sig==True).sum())
    q=bh(np.clip(2*t.p.values,1e-12,1))
    up=((t.p<.025)&(t.score>0)).sum(); dn=((t.p<.025)&(t.score<0)).sum()
    print(f'sheet {s:11} header-cell={hdr:11} rows={len(t)} p<.025={n} IsSig={ns} up={up} down={dn} BH(2p)<.05={(q<.05).sum()} <.10={(q<.10).sum()}')
# ASNS
for s,t in M.items():
    r=t[t.id==76].iloc[0]; print(s,'X76',round(r.score,3),r.p)
# isoprenoid
for tid in (138,140,141):
    for s,t in M.items():
        r=t[t.id==tid].iloc[0]; print('X%d'%tid, r.desc[:35], s, round(r.score,2), r.p, r.sig)
import pickle; pickle.dump(M,open('.arena/run-20261013-005902-s263431/scratch/a008/M4.pkl','wb'))
# fraction down among significant
for s,t in M.items():
    sg=t[t.p<.025]; print(s,'sig',len(sg),'frac down',round((sg.score<0).mean(),2))
# score ranges of strong synergy tasks
