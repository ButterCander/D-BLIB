import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
x = pd.ExcelFile('si/MOESM4.xlsx')
out={}
for s in x.sheet_names:
    d = pd.read_excel(x, sheet_name=s, header=None)
    hdr = d.iloc[0].tolist()
    t = d.iloc[2:, :4].copy(); t.columns=['id','desc','score','p']
    t['id']=pd.to_numeric(t.id,errors='coerce'); t=t.dropna(subset=['id'])
    t['score']=pd.to_numeric(t.score); t['p']=pd.to_numeric(t.p)
    out[s]=t
    print(s, 'hdr', hdr[:3], 'rows', len(t), 'p<.025:', (t.p<0.025).sum(), 'up', ((t.p<0.025)&(t.score>0)).sum(), 'down', ((t.p<0.025)&(t.score<0)).sum(),
          'BH<.05(2p)', (bh(np.clip(2*t.p.values,1e-12,1))<.05).sum(), 'min p', t.p.min())
import pickle; pickle.dump(out, open('.arena/run-20261013-005902-s263431/scratch/a008/m4.pkl','wb'))
# isoprenoid tasks X138,X140,X141 ; ASNS X76
for tid in (138,140,141,76,113,114):
    for s in out:
        r = out[s][out[s].id==tid]
        if len(r): print(tid, s, r.iloc[0].desc[:40], round(r.iloc[0].score,3), r.iloc[0].p)
        else: print(tid, s, 'absent')
