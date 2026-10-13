import openpyxl, pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
B='D:/Cander/School/3-2/system_bio/mtea/'
def sheets(f):
    wb=openpyxl.load_workbook(B+'si/'+f+'.xlsx',read_only=True,data_only=True)
    return {ws.title:list(ws.iter_rows(values_only=True)) for ws in wb.worksheets}
m4=sheets('MOESM4'); m5=sheets('MOESM5'); m7=sheets('MOESM7')
sc={}
for k,rows in m4.items():
    d=pd.DataFrame([r[:5] for r in rows[2:] if r[0] is not None],columns=['id','desc','score','p','sig'])
    d['id']=d.id.astype(int); sc[k]=d.set_index('id')
# BH on authors p (two-sided 2p)
for k in sc:
    p=sc[k].p.astype(float).values
    q=bh(np.clip(2*p,1e-12,1))
    q1=bh(np.clip(p,1e-12,1))
    print(k,'p<.025:',(p<0.025).sum(),'BH(2p)<.05:',(q<.05).sum(),'<.10:',(q<.10).sum(),'| BH(p)<.05:',(q1<.05).sum())
# synergy lists
print()
for k in ['PI3Ki-TAKi','PI3Ki-MEKi']:
    rows=m5[k]; hdr=rows[0]; print(k,hdr[:6],len(rows)-1)
    d=pd.DataFrame([r[:6] for r in rows[1:] if r[0] is not None],columns=['task','desc','A','B','AB','dir'])
    d['id']=d.task.str[1:].astype(int)
    d['strong']=(d.A==0)&(d.B==0)
    print(' rows',len(d),'strong',d.strong.sum(),'weak',(~d.strong).sum(),'up',(d.dir=='up-reg').sum(),'down',(d.dir=='down-reg').sum())
    # margin over best single
    def margin(r):
        if r.dir=='down-reg': return min(r.A,r.B)-r.AB   # positive = more extreme
        return r.AB-max(r.A,r.B)
    d['margin']=d.apply(margin,axis=1)
    w=d[~d.strong]
    print(' weak margin median %.3f  mean %.3f  frac<0.10 %.2f  frac<0.05 %.2f'%(w.margin.median(),w.margin.mean(),(w.margin<0.10).mean(),(w.margin<0.05).mean()))
    # summary sheet
print([r[:3] for r in m5['Summary']])
