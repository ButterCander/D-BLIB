import pandas as pd, numpy as np
exec(open('.arena/run-20261013-005902-s263431/scratch/a001/m4b.py').read().split("for sh in")[0])
for sh in ['TAKi','MEKi','PI3Ki','PI3Ki-MEKi','PI3Ki-TAKi']:
    a=m4(sh); b=m7(sh)
    col=[c for c in b.columns if c.startswith('tide_') and not c.startswith('tide_e')][0]
    j=a.join(b[[col]],how='inner'); j[col]=pd.to_numeric(j[col],errors='coerce')
    sig=j.p<0.025
    print(sh,'n m7',len(b),'m7 nonzero',(j[col]!=0).sum(),'sig4',sig.sum(),'m7 nonzero&sig4',((j[col]!=0)&sig).sum())
    z=j[j[col]!=0]
    print('   corr on m7-nonzero',np.corrcoef(z.score,z[col])[0,1],'exact equal count',(np.isclose(z.score,z[col])).sum(),'of',len(z))
    # among the sig in m4
    s=j[sig]
    print('   among m4-sig: corr',np.corrcoef(s.score,s[col])[0,1],'exact',np.isclose(s.score,s[col]).sum(),'of',len(s))
