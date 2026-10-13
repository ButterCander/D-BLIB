import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
S='D:/Cander/School/3-2/system_bio/mtea/si/'
sheets=['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']
A={}
for s in sheets:
    raw=pd.read_excel(S+'MOESM4.xlsx',sheet_name=s,header=None)
    hdr=raw.iloc[0,2]
    d=raw.iloc[2:,:5].copy(); d.columns=['id','desc','score','p','sig']
    d=d[pd.to_numeric(d.id,errors='coerce').notna()].copy()
    d['id']=d.id.astype(float).astype(int); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    A[s]=d
    n=len(d); sig=(d.p<0.025).sum(); flagged=(d.sig==True).sum()
    q=bh(np.clip(2*d.p.values,1e-12,1))
    print(f'{s:11} header_cell={hdr:11} rows={n} p<.025={sig} flagTrue={flagged} up={((d.p<.025)&(d.score>0)).sum()} down={((d.p<.025)&(d.score<0)).sum()} down%={((d.p<.025)&(d.score<0)).sum()/sig:.3f} BH(2p)<.05={(q<.05).sum()} <.10={(q<.10).sum()} minp={d.p.min()}')
tot=sum((A[s].p<.025).sum() for s in sheets); dn=sum(((A[s].p<.025)&(A[s].score<0)).sum() for s in sheets)
print('total sig',tot,'down',dn,dn/tot)
# text counts: TAKi30 MEKi17 PI77 PIMEK66 PITAK55 -> sum 245
print('paper-text counts sum',30+17+77+66+55)
import pickle; pickle.dump(A,open('D:/Cander/School/3-2/system_bio/mtea/.arena/run-20261013-005902-s263431/scratch/a002/w/A.pkl','wb'))
# tasks of interest
for t in (76,138,140,141,142,113,114,8,45,21,7):
    r=[]
    for s in ['PI3Ki','MEKi','TAKi','PI3Ki-TAKi','PI3Ki-MEKi']:
        x=A[s][A[s].id==t]
        r.append(f'{s}:{x.score.values[0]:+.2f}(p={x.p.values[0]:.4f})' if len(x) else f'{s}:NA')
    desc=A['PI3Ki'][A['PI3Ki'].id==t].desc.values
    print(t,desc[0][:45] if len(desc) else '', ' | '.join(r))
