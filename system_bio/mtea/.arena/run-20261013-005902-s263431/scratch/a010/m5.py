import sys, pickle, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
b=base+'.arena/run-20261013-005902-s263431/scratch/a010/'
x=pd.ExcelFile(base+'si/MOESM5.xlsx')
print(pd.read_excel(x,'Summary',header=None).to_string())
d=pickle.load(open(b+'m4.pkl','rb'))
for sh in ['PI3Ki-TAKi','PI3Ki-MEKi']:
    df=pd.read_excel(x,sh)
    df=df.iloc[:, :6]
    df.columns=['task_id','desc','A','B','AB','dir']
    df=df.dropna(subset=['task_id'])
    print(sh,len(df))
    print(df.dir.value_counts().to_dict())
    other='TAKi' if 'TAKi' in sh else 'MEKi'
    df['strong']=(df.A==0)&(df.B==0)
    # weak margin: |AB| - max(|A|,|B|) in direction
    df['margin']=np.where(df.AB<0, np.minimum(df.A,df.B)-df.AB, df.AB-np.maximum(df.A,df.B))
    w=df[~df.strong]
    print(' strong',df.strong.sum(),'weak',len(w),' weak margin median',w.margin.median(),' frac<0.1',(w.margin<0.1).mean(), 'frac <0.05',(w.margin<0.05).mean())
    print(' all-calls margin median',df.margin.median())
    print(df[['task_id','desc','A','B','AB','dir','margin']].sort_values('margin').head(5).to_string())
    print(' AB score range of strong:',df[df.strong].AB.min(),df[df.strong].AB.max())
    print(df[df.strong][['task_id','desc','AB','dir']].to_string())
    # significance of single drugs in MOESM4: are A,B nonzero values for p-significance?
