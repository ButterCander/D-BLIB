import pandas as pd, numpy as np
B='D:/Cander/School/3-2/system_bio/mtea/'
conds={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki+TAKi':'PI103+5Z','PI3Ki+MEKi':'PI103+PD'}
de={k:pd.read_csv(B+'out/de_%s.tsv'%v,sep='\t',index_col=0) for k,v in conds.items()}
genes=['HMGCS1','HMGCR','MVD','MVK','SCD','SREBF1','FDPS','IDI1','INSIG1','LDLR','ASNS','SREBF2','SCAP']
print('shrunk lfc (padj<.05 marked *)  | unshrunk')
for g in genes:
    s=[];u=[]
    for k in conds:
        r=de[k].loc[g]
        s.append('%+.2f%s'%(r.lfc_shrunk,'*' if r.padj<0.05 else ''))
        u.append('%+.2f'%r.log2FoldChange)
    print('%-7s'%g,' '.join(s),' | ',' '.join(u))
# DEG counts |lfc|>=1 padj<.05 (shrunk and unshrunk)
for k in conds:
    r=de[k]
    print(k,'DEG unshrunk',((r.padj<.05)&(r.log2FoldChange.abs()>=1)).sum(),'shrunk',((r.padj<.05)&(r.lfc_shrunk.abs()>=1)).sum(),'genes',len(r),'padj notNA',r.padj.notna().sum())
