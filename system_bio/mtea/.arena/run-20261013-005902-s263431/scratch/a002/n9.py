import pandas as pd, numpy as np
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv'); gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
print(off.columns.tolist()); print(off.Contrast.unique()[:60])
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SQLE','INSIG1','SREBF2','SCAP','LDLR','SCD','SREBF1']
for tp in ['1h','2h','4h','8h','24h']:
    print('---',tp)
    rows=[]
    for c in ['PI','PD','OXO','PIPD','PIOXO']:
        for pre in ['Diff_','']:
            nm=f'{pre}{c}_{tp}'
            if nm in set(off.Contrast): break
        x=off[off.Contrast==nm].drop_duplicates('sym').set_index('sym')
        rows.append(pd.Series({g:x.logFC.get(g,np.nan) for g in genes},name=c+'('+nm+')'))
    print(pd.DataFrame(rows).round(2).to_string())
