import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
off=pd.read_csv('tsirvouli_official/DEGs_perContrast_all.csv'); gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
off['sym']=gi.reindex(off.Gene_name)['gene_name'].values
off=off.dropna(subset=['sym'])
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SQLE','LDLR','INSIG1','SREBF2','SCAP','RRM1','RRM2','CTPS1']
tp=['1h','2h','4h','8h','24h']
cond=['PI','PD','OXO','PIPD','PIOXO']
for g in ['HMGCS1','HMGCR','MVK','MVD','SQLE','LDLR','SCD','FDPS']:
    print('==',g)
    rows=[]
    for c in cond:
        r=[c]
        for t in tp:
            x=off[(off.Contrast==f'Diff_{c}_{t}')&(off.sym==g)]
            if len(x): 
                x=x.iloc[0]; r.append(f"{x.logFC:+.2f}{'*' if x['adj.P.Val']<0.05 else ''}")
            else: r.append('NA')
        rows.append(r)
    print(pd.DataFrame(rows,columns=['cond']+tp).to_string(index=False))
# program mean 24 genes per cond at 24h and 8h
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
for t in tp:
    print(t,{c: round(off[(off.Contrast==f'Diff_{c}_{t}')&(off.sym.isin(SREBP))].drop_duplicates('sym').logFC.mean(),2) for c in cond})
