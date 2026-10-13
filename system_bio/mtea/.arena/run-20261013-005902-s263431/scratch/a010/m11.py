import sys, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
off=pd.read_csv(base+'tsirvouli_official/DEGs_perContrast_all.csv')
print(off.columns.tolist(), off.shape)
print(off.Contrast.unique()[:60])
gi=pd.read_csv(base+'tsirvouli/gene_info.tsv',sep='\t',index_col=0)
print(gi.head(3))
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
genes=['HMGCS1','HMGCR','MVD','MVK','IDI1','FDPS','SCD','SQLE','LDLR','INSIG1','SREBF2','SCAP','ASNS']
for tp in ['1h','2h','4h','8h','24h']:
    print('time',tp)
    rows=[]
    for g in genes:
        r=[g]
        for c in ['PI','PD','OXO','PIPD','PIOXO']:
            x=off[(off.Contrast==f'Diff_{c}_{tp}')&(off.sym==g)]
            if len(x)==0: r.append('NA')
            else:
                x=x.iloc[0]
                pcol=[k for k in off.columns if k.lower() in ('adj.p.val','adj_p_val','padj','adj.p','fdr')]
                r.append(f"{x.logFC:+.2f}" + ('*' if pcol and x[pcol[0]]<0.05 else ''))
        rows.append(r)
    print(pd.DataFrame(rows,columns=['gene','PI','PD','5Z','PI+PD','PI+5Z']).to_string(index=False))
