import sys, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
LINES=["CAPAN2","CFPAC1","MiaPaca2","PANC1"]
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SQLE','LDLR','INSIG1','RRM1','RRM2','CTPS1']
D={(l,t):pd.read_csv(base+f'out/pdac_de_{l}_{t}.tsv',sep='\t').set_index('gene') for l in LINES for t in ('h3','h12','h24')}
print('24h trametinib logFC (adj.P<.05 *)')
rows=[]
for g in genes:
    r=[g]
    for l in LINES:
        x=D[(l,'h24')].loc[g] if g in D[(l,'h24')].index else None
        r.append('NA' if x is None else f"{x.logFC:+.2f}{'*' if x['adj.P.Val']<0.05 else ''}")
    rows.append(r)
print(pd.DataFrame(rows,columns=['gene']+LINES).to_string(index=False))
print()
for t in ('h3','h12','h24'):
    print(t,{l:round(D[(l,t)]['logFC'].reindex([g for g in SREBP if g in D[(l,t)].index]).mean(),2) for l in LINES})
print('n genes up>=0.5 & padj<.05 among 24 in each line at 24h:')
for l in LINES:
    d=D[(l,'h24')].reindex(SREBP).dropna()
    print(l, 'n',len(d),'up>0',(d.logFC>0).sum(),'>=0.5 & sig',((d.logFC>=0.5)&(d['adj.P.Val']<0.05)).sum())
# HMGCR >=0.5 in how many lines
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD']:
    n=sum(1 for l in LINES if (D[(l,'h24')].loc[g].logFC>=0.5 and D[(l,'h24')].loc[g]['adj.P.Val']<0.05))
    print(g,'lines with logFC>=0.5 & padj<.05 at 24h:',n)
