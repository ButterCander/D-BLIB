import openpyxl, pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
B='D:/Cander/School/3-2/system_bio/mtea/'
def sheets(f):
    wb=openpyxl.load_workbook(B+'si/'+f+'.xlsx',read_only=True,data_only=True)
    return {ws.title:list(ws.iter_rows(values_only=True)) for ws in wb.worksheets}
m4=sheets('MOESM4'); m5=sheets('MOESM5')
sc={}
for k,rows in m4.items():
    d=pd.DataFrame([r[:5] for r in rows[2:] if r[0] is not None],columns=['id','desc','score','p','sig'])
    d['id']=d.id.astype(int); sc[k]=d.set_index('id')
for comb,single in [('PI3Ki-TAKi','TAKi'),('PI3Ki-MEKi','MEKi')]:
    rows=m5[comb]
    d=pd.DataFrame([r[:6] for r in rows[1:] if r[0] is not None],columns=['task','desc','A','B','AB','dir'])
    d['id']=d.task.str[1:].astype(int)
    st=d[(d.A==0)&(d.B==0)]
    print(comb,'strong calls:')
    for _,r in st.iterrows():
        pA=sc['PI3Ki'].p[r.id]; pB=sc[single].p[r.id]; pAB=sc[comb].p[r.id]
        print('  X%d %-45s AB=%.3f pPI3Ki=%.4f p%s=%.4f pAB=%.4f'%(r.id,r.desc[:45],r.AB,pA,single,pB,pAB))
    # all calls: AB vs control q
    q=pd.Series(bh(np.clip(2*sc[comb].p.astype(float).values,1e-12,1)),index=sc[comb].index)
    ids=d.id.values
    print('  calls:',len(ids),' q<.05:',(q[ids]<.05).sum(),' q<.10:',(q[ids]<.10).sum())
    print('  strong q:',[round(q[i],3) for i in st.id])
# Table: X138,X140,X141, MEKi sheet
for t in [138,140,141]:
    print('X%d'%t, {k:(round(float(sc[k].score[t]),2),sc[k].p[t]) for k in sc}, sc['MEKi'].desc[t])
