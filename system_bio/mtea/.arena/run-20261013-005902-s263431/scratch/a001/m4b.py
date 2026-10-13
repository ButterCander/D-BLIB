import pandas as pd, numpy as np
def m4(sh):
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=sh,header=None).iloc[2:,:4]; d.columns=['id','desc','score','p']
    d=d.dropna(subset=['id']); d['id']='X'+d.id.astype(int).astype(str); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    return d.set_index('id')
def m7(sh):
    d=pd.read_excel('si/MOESM7.xlsx',sheet_name=sh)
    d=d.set_index('task_id')
    return d
for sh in ['TAKi','MEKi','PI3Ki','PI3Ki-MEKi','PI3Ki-TAKi']:
    a=m4(sh); b=m7(sh)
    col=[c for c in b.columns if c.startswith('tide_') and not c.startswith('tide_e')][0]
    j=a.join(b[[col]],how='inner')
    j[col]=pd.to_numeric(j[col],errors='coerce')
    print(sh,len(j),'corr sheet score vs MOESM7',col,np.corrcoef(j.score,j[col])[0,1], 'max abs diff',(j.score-j[col]).abs().max())
# cross: TAKi sheet vs MOESM7 MEKi
a=m4('TAKi'); b=m7('MEKi'); col=[c for c in b.columns if c.startswith('tide_') and not c.startswith('tide_e')][0]
j=a.join(b[[col]],how='inner'); j[col]=pd.to_numeric(j[col],errors='coerce'); print('TAKi sheet vs MOESM7 MEKi',np.corrcoef(j.score,j[col])[0,1])
