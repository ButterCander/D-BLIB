import pandas as pd, numpy as np
xl = pd.ExcelFile('si/MOESM4.xlsx')
S = {}
for sh in xl.sheet_names:
    d = pd.read_excel(xl, sh, header=None)
    print(sh, d.shape, d.iloc[0,:3].tolist())
    t = d.iloc[2:, :4].copy(); t.columns=['id','desc','score','p']
    t = t.dropna(subset=['id'])
    t['id']=t.id.astype(float).astype(int); t['score']=pd.to_numeric(t.score); t['p']=pd.to_numeric(t.p)
    S[sh]=t
    print('  n rows', len(t), 'p<0.025:', (t.p<0.025).sum(), ' p<=0.025:', (t.p<=0.025).sum(), 'neg among sig:', ((t.p<0.025)&(t.score<0)).sum(), 'pos', ((t.p<0.025)&(t.score>0)).sum(), 'zero scores', (t.score==0).sum(), 'min/max p', t.p.min(), t.p.max())
import pickle; pickle.dump(S, open('.arena/run-20261013-005902-s263431/scratch/a004/moesm4.pkl','wb'))
