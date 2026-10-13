import pandas as pd, numpy as np, pickle
A=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a002/A.pkl','rb'))
for t in (138,140,141,142,76,113,112,178,8,17,68):
    row=[]
    for k in ['PI3Ki','MEKi','TAKi','PI3Ki-TAKi','PI3Ki-MEKi']:
        d=A[k]; r=d[d.id==t]
        row.append(f"{k}: {r.score.values[0]:+.2f} (p={r.p.values[0]:.4f})")
    print(t, A['PI3Ki'][A['PI3Ki'].id==t].desc.values[0][:40],'|',' | '.join(row))
