import pandas as pd, numpy as np
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv')
print(off.columns.tolist(), off.shape)
print(off.Contrast.unique()[:60])
gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
print(gi.columns.tolist()); print(gi.head(3))
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
print(off.sym.notna().mean())
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','INSIG1','LDLR','SCD','SREBF2','SCAP','SQLE']
for t in ['1h','4h','8h','24h']:
    print('--',t)
    for g in genes:
        r=[]
        for c in ['PI','PD','OXO','PIPD','PIOXO']:
            x=off[(off.Contrast==f'Diff_{c}_{t}')&(off.sym==g)]
            r.append('%s %+.2f%s'%(c,x.logFC.iloc[0],'*' if x['adj.P.Val'].iloc[0]<0.05 else '') if len(x) else c+' NA')
        print('%-7s'%g,' | '.join(r))
