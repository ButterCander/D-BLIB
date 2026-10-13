import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv')
gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
print(off.Gene_name.head(3).tolist())
off['sym']=off.Gene_name.map(gi['gene_name']) if off.Gene_name.iloc[0].startswith('ENSG') else off.Gene_name
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','SCD','LDLR','INSIG1','SREBF2','SCAP','SREBF1']
for t in ['1h','4h','8h','24h']:
    print('--',t)
    for g in genes:
        out=[]
        for c in ['PI','PD','OXO','PIPD','PIOXO']:
            r=off[(off.Contrast==f'Diff_{c}_{t}')&(off.sym==g)]
            out.append(f'{c}:{r.logFC.iloc[0]:+.2f}({r["adj.P.Val"].iloc[0]:.3f})' if len(r) else f'{c}:NA')
        print(f'  {g:7s}',' '.join(out))
