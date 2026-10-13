import pandas as pd, numpy as np
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
x=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
for p,name in [('A','AsPC-1'),('P','PANC-1'),('C','CFPAC-1'),('S9','SW1990')]:
    lf=np.log2(x[[f'{p}_M_{i}' for i in (1,2,3)]]+1).mean(axis=1)-np.log2(x[[f'{p}_NT_{i}' for i in (1,2,3)]]+1).mean(axis=1)
    gs=[g for g in SREBP if g in lf.index]
    print(name,'program mean %.2f'%lf[gs].mean(),'n genes',len(gs),'n up',int((lf[gs]>0).sum()),{g:round(lf[g],2) for g in ['HMGCS1','HMGCR','MVD','MVK','SCD'] if g in lf.index}, 'total genes in table',len(x))
# statin gene-level
for cl in ['PANC1','MiaPaca2']:
    for st in ['Simvastatin','Fluvastatin','Lovastatin']:
        d=pd.read_csv(f'out/statin_de_{cl}_{st}_High.tsv',sep='\t').set_index('gene')
        print(cl,st,{g:(round(d.loc[g,'logFC'],2)) for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','LDLR','SQLE'] if g in d.index})
# 24-gene program in statin
for cl in ['PANC1','MiaPaca2']:
    d=pd.read_csv(f'out/statin_de_{cl}_Simvastatin_High.tsv',sep='\t').set_index('gene')
    v=d.reindex([g for g in SREBP if g in d.index]).logFC; print(cl,'simva high program',round(v.mean(),2),'n up',int((v>0).sum()),'of',len(v))
