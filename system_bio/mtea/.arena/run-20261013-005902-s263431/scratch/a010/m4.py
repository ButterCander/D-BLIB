import sys, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
x=pd.ExcelFile(base+'si/MOESM4.xlsx')
d={}
for sh in x.sheet_names:
    df=pd.read_excel(x,sh,header=1)
    hdr=pd.read_excel(x,sh,header=None,nrows=1).iloc[0].tolist()
    df=df.rename(columns={'Task ID':'id','TIDE score':'score','p-value':'p','Is Significant':'sig'})
    d[sh]=df[['id','Description','score','p','sig']].dropna(subset=['id'])
    print(sh,'header cell',hdr,'n rows',len(d[sh]),'n p<0.025',(d[sh].p<0.025).sum(),'n sig True',d[sh].sig.sum(),'down',((d[sh].p<0.025)&(d[sh].score<0)).sum(),'up',((d[sh].p<0.025)&(d[sh].score>0)).sum())
import pickle
pickle.dump(d,open(base+'.arena/run-20261013-005902-s263431/scratch/a010/m4.pkl','wb'))
# MOESM7 compare
x7=pd.ExcelFile(base+'si/MOESM7.xlsx')
print(x7.sheet_names)
s=pd.read_excel(x7,'Summary',header=None)
print(s.head(12).to_string())
