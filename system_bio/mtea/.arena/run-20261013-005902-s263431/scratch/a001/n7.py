import pandas as pd, numpy as np
conds=[('PI3Ki','PI103'),('MEKi','PD0325901'),('TAKi','5Z'),('PI3Ki+TAKi','PI103+5Z'),('PI3Ki+MEKi','PI103+PD')]
D={n:pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0) for n,f in conds}
print(D['PI3Ki'].columns.tolist())
for g in ['HMGCS1','HMGCR','MVD','MVK','IDI1','FDPS','SCD','SREBF1','SREBF2','INSIG1','LDLR','SQLE','LSS','ASNS']:
    r=[]
    for n,_ in conds:
        d=D[n]
        lf=d.loc[g,'lfc_shrunk']; lu=d.loc[g,'log2FoldChange']; p=d.loc[g,'padj']
        r.append('%+.2f%s(%+.2f)'%(lf,'*' if p<0.05 else '',lu))
    print('%-7s'%g,' | '.join(r))
