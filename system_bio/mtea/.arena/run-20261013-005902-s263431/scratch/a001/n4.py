import pandas as pd, numpy as np
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
lines=['CAPAN2','CFPAC1','MiaPaca2','PANC1']
rows={}
for l in lines:
    for h in ['h3','h12','h24']:
        d=pd.read_csv(f'out/pdac_de_{l}_{h}.tsv',sep='\t').set_index('gene')
        rows[(l,h)]=d
print('program mean (24 listed genes) per line/time')
for l in lines:
    out=[]
    for h in ['h3','h12','h24']:
        d=rows[(l,h)]['logFC']; g=[s for s in SREBP if s in d.index]
        out.append('%s n=%d %.2f'%(h,len(g),d.reindex(g).mean()))
    print(l,out)
print()
genes=['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SREBF1','SREBF2','INSIG1','LDLR','SQLE','LSS','DHCR7']
t=[]
for g in genes:
    r={}
    for l in lines:
        d=rows[(l,'h24')]
        if g in d.index: r[l]='%+.2f%s'%(d.loc[g,'logFC'],'*' if d.loc[g,'adj.P.Val']<0.05 else '')
        else: r[l]='NA'
    t.append((g,r))
for g,r in t: print(g,r)
print({l:int((rows[(l,'h24')]['adj.P.Val']<0.05).sum()) for l in lines},'DEG padj<.05 at 24h')
