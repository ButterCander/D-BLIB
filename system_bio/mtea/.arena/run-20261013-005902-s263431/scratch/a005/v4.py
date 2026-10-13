import pandas as pd, numpy as np, json
from scipy import stats
pd.set_option('display.width',250)
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
print(len(SREBP))
res={}
for l in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    out=[]
    for h in ['h3','h12','h24']:
        d=pd.read_csv(f'out/pdac_de_{l}_{h}.tsv',sep='\t').set_index('gene')
        pm=d.reindex([g for g in SREBP if g in d.index]).logFC.mean()
        out.append((h,round(pm,2),int((d['adj.P.Val']<0.05).sum()),len([g for g in SREBP if g in d.index])))
    print(l,out)
# 4 trunk genes in each line 24h with padj
for l in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    d=pd.read_csv(f'out/pdac_de_{l}_h24.tsv',sep='\t').set_index('gene')
    print(l,{g:(round(d.loc[g,'logFC'],2),'*' if d.loc[g,'adj.P.Val']<.05 else 'ns') for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','IDI1','FDPS','SQLE','LDLR','INSIG1'] if g in d.index})
# claim: HMGCS1, MVK, MVD, SCD up >=0.5 in 3-4 of 4 lines
cnt={}
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD']:
    n=0
    for l in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
        d=pd.read_csv(f'out/pdac_de_{l}_h24.tsv',sep='\t').set_index('gene')
        if g in d.index and d.loc[g,'logFC']>=0.5 and d.loc[g,'adj.P.Val']<0.05: n+=1
    cnt[g]=n
print('lines with logFC>=0.5 & padj<.05 (of 4):',cnt)
# E-flux
g=json.load(open('out/eflux_grid.json'))
vals=[v for x in g for v in x['pred'].values()]
print('eflux settings',len(g),'min',min(vals),'max',max(vals), 'n nonzero settings', sum(any(abs(v)>1e-6 for v in x['pred'].values()) for x in g))
print([ (x['vmax'],x['floor'],{k:round(v,1) for k,v in x['pred'].items()},round(x['rho'],2)) for x in g][:18])
