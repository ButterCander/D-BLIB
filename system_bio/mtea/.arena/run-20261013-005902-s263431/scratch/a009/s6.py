import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TRUNK = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
LINES=["CAPAN2","CFPAC1","MiaPaca2","PANC1"]
D={(l,t):pd.read_csv(f'out/pdac_de_{l}_{t}.tsv',sep='\t').set_index('gene') for l in LINES for t in ('h3','h12','h24')}
print('--- 24h gene-level logFC (padj)')
for gn in ['HMGCS1','HMGCR','MVK','MVD','SCD','IDI1','FDPS','SQLE','LDLR','INSIG1']:
    s=[]
    for l in LINES:
        d=D[(l,'h24')]
        if gn in d.index: s.append(f"{l}:{d.loc[gn,'logFC']:+.2f}({d.loc[gn,'adj.P.Val']:.3f})")
        else: s.append(f'{l}:NA')
    print(gn,' '.join(s))
print('--- program means (24 genes) and trunk(6)')
for l in LINES:
    out=[]
    for t in ('h3','h12','h24'):
        d=D[(l,t)]['logFC']
        out.append(f"{t}: prog {d.reindex([g for g in SREBP if g in d.index]).mean():+.2f} trunk {d.reindex([g for g in TRUNK if g in d.index]).mean():+.2f}")
    print(l,' | '.join(out))
print('--- DEG counts padj<0.05 at 24h (all |lfc|)')
for l in LINES:
    d=D[(l,'h24')]; print(l,(d['adj.P.Val']<0.05).sum(), 'n genes',len(d))
# 4-gene (HMGCS1,HMGCR,MVK,MVD) mean at 24h
for l in LINES:
    d=D[(l,'h24')]['logFC']; print(l,'4gene mean',d.reindex(['HMGCS1','HMGCR','MVK','MVD']).mean().round(2))
# count lines with logFC>=0.5 & padj<.05 per gene
print('--- #lines (of 4) with lfc>=0.5 & padj<0.05 at 24h')
for gn in ['HMGCS1','HMGCR','MVK','MVD','SCD','IDI1','FDPS']:
    n=sum(1 for l in LINES if gn in D[(l,'h24')].index and D[(l,'h24')].loc[gn,'logFC']>=0.5 and D[(l,'h24')].loc[gn,'adj.P.Val']<0.05)
    n2=sum(1 for l in LINES if gn in D[(l,'h24')].index and D[(l,'h24')].loc[gn,'logFC']>0 and D[(l,'h24')].loc[gn,'adj.P.Val']<0.05)
    print(gn,n,'>=0.5 ;',n2,'any up sig')
