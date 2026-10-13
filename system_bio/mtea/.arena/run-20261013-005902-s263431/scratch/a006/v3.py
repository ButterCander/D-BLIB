import pandas as pd, numpy as np, openpyxl
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=True)
S={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    df=pd.DataFrame([r[:4] for r in rows[2:]],columns=['id','desc','score','p']).dropna(subset=['id'])
    df['task']='X'+df.id.astype(int).astype(str); df['score']=pd.to_numeric(df.score); df['p']=pd.to_numeric(df.p)
    S[ws.title]=df.set_index('task')
for t in ['X138','X139','X140','X141','X76','X8','X124']:
    print(t, S['PI3Ki'].loc[t,'desc'])
    for s in ['PI3Ki','MEKi','TAKi','PI3Ki-TAKi','PI3Ki-MEKi']:
        r=S[s].loc[t]; print('   ',s, round(r.score,3), r.p)
# ASNS task: PI3Ki -5.9 ; rank
for s in S:
    d=S[s]; print(s,'X76 score',round(d.loc['X76','score'],2),'p',d.loc['X76','p'],'rank of most negative',int(d.score.rank().loc['X76']))
# up/down counts among significant
for s in S:
    d=S[s]; sig=d[d.p<0.025]; print(s,'sig',len(sig),'up',(sig.score>0).sum(),'down',(sig.score<0).sum(), 'zero', (sig.score==0).sum())
