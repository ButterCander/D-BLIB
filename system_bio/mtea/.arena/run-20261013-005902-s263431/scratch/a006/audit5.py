import pandas as pd, numpy as np
B='D:/Cander/School/3-2/system_bio/mtea/'
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
lines=['CAPAN2','CFPAC1','MiaPaca2','PANC1']
print('24 h trametinib (GSE238198) per-gene logFC (adj.P<.05 *)')
tab={}
for l in lines:
    d=pd.read_csv(B+'out/pdac_de_%s_h24.tsv'%l,sep='\t').set_index('gene')
    tab[l]=d
    prog=d.reindex([s for s in SREBP if s in d.index])
    print(l,'program mean %.2f  n=%d  up=%d  sig up=%d  | DEGs padj<.05: %d'%(prog.logFC.mean(),prog.logFC.notna().sum(),(prog.logFC>0).sum(),((prog.logFC>0)&(prog['adj.P.Val']<.05)).sum(),(d['adj.P.Val']<.05).sum()))
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','IDI1','FDPS','SQLE','LDLR','INSIG1','RRM1','RRM2','CTPS1']:
    print('%-7s'%g,' '.join('%s %+.2f%s'%(l,tab[l].loc[g,'logFC'],'*' if tab[l].loc[g,'adj.P.Val']<.05 else '') for l in lines))
# proportion of lines >=0.5
for g in ['HMGCS1','HMGCR','MVK','MVD','SCD']:
    print(g,'lines logFC>=0.5:',sum(tab[l].loc[g,'logFC']>=0.5 for l in lines),' sig&>=0.5:',sum((tab[l].loc[g,'logFC']>=0.5) and (tab[l].loc[g,'adj.P.Val']<.05) for l in lines))
# time course
for l in lines:
    out=[]
    for h in ['h3','h12','h24']:
        d=pd.read_csv(B+'out/pdac_de_%s_%s.tsv'%(l,h),sep='\t').set_index('gene')
        out.append('%s %+.2f'%(h,d.reindex([s for s in SREBP if s in d.index]).logFC.mean()))
    print(l,out)
