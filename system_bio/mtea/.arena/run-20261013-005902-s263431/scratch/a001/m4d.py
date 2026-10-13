import pandas as pd, numpy as np
exec(open('.arena/run-20261013-005902-s263431/scratch/a001/m4b.py').read().split("for sh in")[0])
for sh in ['TAKi','MEKi','PI3Ki','PI3Ki-MEKi','PI3Ki-TAKi']:
    a=m4(sh); b=m7(sh)
    col=[c for c in b.columns if c.startswith('tide_') and not c.startswith('tide_e')][0]
    j=a.join(b[[col]],how='inner'); j[col]=pd.to_numeric(j[col],errors='coerce')
    sig=j.p<0.025; nz=j[col]!=0
    d1=j[sig&~nz]; d2=j[~sig&nz]
    print(sh,'sig in M4 but 0 in M7:',[(i,r.p) for i,r in d1.iterrows()],' | nonzero in M7 but p>=.025 in M4:',[(i,r.p) for i,r in d2.iterrows()])
