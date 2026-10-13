import sys; sys.path.insert(0,'.arena/run-20261013-005902-s263431/scratch/a004')
from load4 import *
import pandas as pd, numpy as np
d=load4()
mine={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
tm={}
for k,v in mine.items():
    t=pd.read_csv(f'out/tide_p10000_{v}.tsv',sep='\t'); t['tid']=t.task_id.str.replace('X','').astype(int); tm[k]=t.set_index('tid')
for sheet,(hdr,df) in d.items():
    df=df.set_index('tid')
    row=[]
    for k in mine:
        j=df.join(tm[k],how='inner',lsuffix='_x',rsuffix='_y')
        row.append((k,round(np.corrcoef(j.score_x,j.score_y)[0,1],3),round(j[['score_x','score_y']].corr(method='spearman').iloc[0,1],3)))
    print(sheet,hdr,row)
# BH on 2p
for sheet,(hdr,df) in d.items():
    df=df.dropna(subset=['p']); q=bh(np.minimum(2*df.p.values,1))
    print(sheet,'p<.025',(df.p<0.025).sum(),'BH(2p)<.05',(q<0.05).sum(),'<.10',(q<0.10).sum())
