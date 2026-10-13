import openpyxl, numpy as np, pandas as pd
from scipy import stats
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=True)
S={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    hdr=rows[0][2]
    d=[(r[0],r[1],r[2],r[3]) for r in rows[2:] if r[0] is not None and isinstance(r[0],(int,float)) and r[3] is not None]
    df=pd.DataFrame(d,columns=['id','desc','score','p'])
    S[ws.title]=df
    print(ws.title,'headercell=',hdr,'n=',len(df))
def bh(p):
    p=np.asarray(p);n=len(p);o=np.argsort(p);q=np.empty(n);r=p[o]*n/(np.arange(n)+1)
    r=np.minimum.accumulate(r[::-1])[::-1];q[o]=np.minimum(r,1);return q
mine={'TAKi':'5Z','MEKi':'PD0325901','PI3Ki':'PI103','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
for k,df in S.items():
    q=bh(df.p.values)
    ndn=(df.p<0.025).sum()
    nsig=df['score'][df.p<0.025]
    print(k,'p<.025:',ndn,'BH<.05:',(q<0.05).sum(),'BH<.10:',(q<0.10).sum(),'neg frac among sig:',round((nsig<0).mean(),2))
    t=pd.read_csv('out/tide_%s.tsv'%mine[k],sep='\t')
    t['num']=t.task_id.str.replace('X','').astype(int)
    m=df.merge(t,left_on='id',right_on='num')
    print('   r(sheet score, my score) =',round(stats.pearsonr(m.score_x,m.score_y)[0],3),' n=',len(m))
# cross: sheet TAKi vs my PD, sheet MEKi vs my 5Z
for k,o in [('TAKi','PD0325901'),('MEKi','5Z')]:
    t=pd.read_csv('out/tide_%s.tsv'%o,sep='\t');t['num']=t.task_id.str.replace('X','').astype(int)
    m=S[k].merge(t,left_on='id',right_on='num');print('cross',k,o,round(stats.pearsonr(m.score_x,m.score_y)[0],3))
import pickle;pickle.dump(S,open('.arena/run-20261013-005902-s263431/scratch/a005/S.pkl','wb'))
for k in S:
    d=S[k]
    for tid in (138,140,141):
        r=d[d.id==tid]
        print(k,tid,r.desc.values[0][:30],round(r.score.values[0],3),r.p.values[0])
