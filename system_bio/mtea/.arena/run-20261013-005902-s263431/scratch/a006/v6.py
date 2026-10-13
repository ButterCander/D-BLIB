import pandas as pd, numpy as np, openpyxl
from scipy.stats import false_discovery_control as bh
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=True)
S={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    df=pd.DataFrame([r[:4] for r in rows[2:]],columns=['id','desc','score','p']).dropna(subset=['id'])
    df['task']='X'+df.id.astype(int).astype(str); df['score']=pd.to_numeric(df.score); df['p']=pd.to_numeric(df.p)
    S[ws.title]=df.set_index('task')
w5=openpyxl.load_workbook('si/MOESM5.xlsx',read_only=True,data_only=True)
for ws in w5.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    if ws.title=='Summary':
        for r in rows: print(r[:4])
        continue
    df=pd.DataFrame([r[:6] for r in rows[1:]],columns=['task','desc','A','B','AB','dir']).dropna(subset=['task'])
    print(ws.title,'rows',len(df),'direction counts',df.dir.value_counts().to_dict())
    other='TAKi' if ws.title=='PI3Ki-TAKi' else 'MEKi'
    # classify strong: both singles zero
    strong=df[(df.A==0)&(df.B==0)]; weak=df[~((df.A==0)&(df.B==0))]
    print('  strong',len(strong),'weak',len(weak))
    # BH q for combination vs control
    p=S[ws.title].p; q=pd.Series(bh(np.clip(2*p.values,1e-12,1)),index=p.index)
    for nm,d in (('strong',strong),):
        for t,r in d.iterrows():
            tid=r.task
            print('   ',tid,r.desc[:45],'AB score',round(r.AB,2),'p_AB',S[ws.title].p.get(tid),'q_AB(2p BH)',round(q.get(tid),3),'p_PI',S['PI3Ki'].p.get(tid),'p_other',S[other].p.get(tid))
    # weak margins
    wk=weak.copy(); wk['margin']=np.where(wk.dir=='down-reg', wk[['A','B']].min(axis=1)-wk.AB, wk.AB-wk[['A','B']].max(axis=1))
    print('  weak margin median',round(wk.margin.median(),3),'frac<0.10',round((wk.margin<0.10).mean(),2),'n',len(wk))
    # how many of listed tasks are significant in AB (p<.025)?
    print('  listed tasks with p_AB<.025:',int((S[ws.title].p.reindex(df.task)<0.025).sum()),'of',len(df))
