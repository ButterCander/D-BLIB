import sys, pickle, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
b=base+'.arena/run-20261013-005902-s263431/scratch/a010/'
d=pickle.load(open(b+'m4.pkl','rb'))
for t in [138,140,141]:
    print('task X%d'%t)
    for k in ['MEKi','TAKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
        r=d[k][d[k].id==t].iloc[0]
        print('   ',k,round(r.score,3),r.p,r.sig, r.Description)
