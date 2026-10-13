import pandas as pd, numpy as np, openpyxl
from scipy.stats import false_discovery_control as bh, spearmanr, pearsonr
pd.set_option('display.width',250)
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=True)
S={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    hdr=rows[0]
    df=pd.DataFrame([r[:4] for r in rows[2:]],columns=['id','desc','score','p'])
    df=df.dropna(subset=['id']); df['p']=pd.to_numeric(df.p,errors='coerce'); df['score']=pd.to_numeric(df.score,errors='coerce')
    S[ws.title]=df
    print(ws.title,'header cell',hdr[2],'rows',len(df),'p<.025:',(df.p<0.025).sum(),'min p',df.p.min(),'max p',df.p.max())
# BH on authors p, 2p
for k,df in S.items():
    p=np.clip(2*df.p.values,1e-12,1)
    q=bh(p)
    print(k,'BH q<.05 (2p):',(q<.05).sum(),' q<.10:',(q<.10).sum(),' tail p BH<.05:',(bh(df.p.values)<.05).sum(), 'n tasks',len(df))
