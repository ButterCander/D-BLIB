import pandas as pd, numpy as np
pd.set_option('display.width',250)
d=pd.read_csv('out/tidex_statin_PANC1_Simvastatin_High.tsv',sep='\t'); print(d.columns.tolist()); print(d.head(3))
for cl in ['PANC1','MiaPaca2']:
    for dose in ['High','Low']:
        for st in ['Fluvastatin','Lovastatin','Simvastatin']:
            t=pd.read_csv(f'out/tidex_statin_{cl}_{st}_{dose}.tsv',sep='\t').set_index('task_id')
            sc=t['score']
            r=sc.rank(ascending=False,method='min')
            n=sc.notna().sum(); npos=(sc>0).sum(); nnz=(sc!=0).sum()
            pcol=[c for c in t.columns if 'p' in c.lower() and 'val' in c.lower()]
            row={k:(int(r[k]),round(float(sc[k]),2), round(float(t.loc[k,pcol[0]]),3) if pcol else None) for k in ['X138','X140','X141']}
            print(cl,dose,st,'n',n,'n>0',npos,'nonzero',nnz,row)
