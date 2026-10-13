import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
L=['CAPAN2','CFPAC1','MiaPaca2','PANC1']
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SQLE','LDLR','INSIG1','SREBF2','RRM1','RRM2','CTPS1','FOSL1']
for h in ['h24','h12','h3']:
    print('==',h)
    rows=[]
    for g in genes:
        r=[g]
        for l in L:
            d=pd.read_csv(f'out/pdac_de_{l}_{h}.tsv',sep='\t').set_index('gene')
            if g in d.index:
                x=d.loc[g]; r.append(f"{x.logFC:+.2f}{'*' if x['adj.P.Val']<0.05 else ''}")
            else: r.append('NA')
        rows.append(r)
    print(pd.DataFrame(rows,columns=['gene']+L).to_string(index=False))
# program mean
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
for h in ['h3','h12','h24']:
    out=[]
    for l in L:
        d=pd.read_csv(f'out/pdac_de_{l}_{h}.tsv',sep='\t').set_index('gene')['logFC']
        s=[g for g in SREBP if g in d.index]
        out.append((l,len(s),round(d.reindex(s).mean(),2),int((d.reindex(s)>0).sum())))
    print(h,out)
