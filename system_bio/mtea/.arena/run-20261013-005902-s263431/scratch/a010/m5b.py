import sys, pickle, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
b=base+'.arena/run-20261013-005902-s263431/scratch/a010/'
d=pickle.load(open(b+'m4.pkl','rb'))
x=pd.ExcelFile(base+'si/MOESM5.xlsx')
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); q=np.empty(n); 
    r=p[o]*n/np.arange(1,n+1); r=np.minimum.accumulate(r[::-1])[::-1]; q[o]=np.minimum(r,1); return q
for sh,single in [('PI3Ki-TAKi','TAKi'),('PI3Ki-MEKi','MEKi')]:
    df=pd.read_excel(x,sh).iloc[:, :6]; df.columns=['id','desc','A','B','AB','dir']
    df['id']=df['id'].str.replace('X','').astype(int)
    comb=d[sh].set_index('id'); pi=d['PI3Ki'].set_index('id'); sg=d[single].set_index('id')
    comb['q2']=bh(np.minimum(1,2*comb.p.values)); 
    df['p_comb']=df.id.map(comb.p); df['q_comb']=df.id.map(comb.q2)
    df['p_pi']=df.id.map(pi.p); df['p_single']=df.id.map(sg.p)
    df['s_pi']=df.id.map(pi.score); df['s_single']=df.id.map(sg.score); df['s_comb']=df.id.map(comb.score)
    df['strong']=(df.A==0)&(df.B==0)
    print(sh)
    print(df[df.strong][['id','desc','AB','p_comb','q_comb','p_pi','p_single','s_pi','s_single']].to_string())
    print('all calls sig in comb sheet:', (df.p_comb<0.025).all(), ' n q<0.05',(df.q_comb<0.05).sum(),' q<.10',(df.q_comb<0.10).sum())
    # check that A,B equal sheet scores where sig else 0
    chk=((df.A-np.where(df.p_pi<0.025,df.s_pi,0)).abs()<1e-9).mean()
    print('A matches PI3Ki sheet masked:',chk)
    chk=((df.B-np.where(df.p_single<0.025,df.s_single,0)).abs()<1e-9).mean()
    print('B matches single sheet masked:',chk)
