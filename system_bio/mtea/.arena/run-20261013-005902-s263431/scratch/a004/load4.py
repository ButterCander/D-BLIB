import openpyxl, pandas as pd, numpy as np
def load4():
    wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=False)
    out={}
    for ws in wb.worksheets:
        rows=list(ws.iter_rows(values_only=True))
        hdr=rows[0]
        recs=[]
        for r in rows[2:]:
            if r[0] is None: continue
            try: tid=int(float(r[0]))
            except: continue
            recs.append((tid,r[1],r[2],r[3]))
        df=pd.DataFrame(recs,columns=['tid','desc','score','p'])
        df['score']=pd.to_numeric(df['score'],errors='coerce'); df['p']=pd.to_numeric(df['p'],errors='coerce')
        out[ws.title]=(hdr[2],df)
    return out
def bh(p):
    p=np.asarray(p,float); n=len(p); o=np.argsort(p); r=np.empty(n); 
    q=p[o]*n/np.arange(1,n+1); q=np.minimum.accumulate(q[::-1])[::-1]; r[o]=np.minimum(q,1); return r
if __name__=='__main__':
    d=load4()
    for k,(h,df) in d.items():
        df=df.dropna(subset=['p'])
        sig=df[df.p<0.025]
        print(k,'hdr cell=',h,'n=',len(df),'sig<.025=',len(sig),'down=',(sig.score<0).sum(),'up=',(sig.score>0).sum(),'BH<.05=',(bh(df.p.values)<0.05).sum(), 'BH<.10=',(bh(df.p.values)<0.10).sum(), 'p==0 n=',(df.p==0).sum())
