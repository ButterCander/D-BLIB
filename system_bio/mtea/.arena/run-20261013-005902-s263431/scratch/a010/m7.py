import sys, pandas as pd, numpy as np
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
conds=[('PI3Ki','PI103'),('MEKi','PD0325901'),('TAKi','5Z'),('PI3Ki+TAKi','PI103+5Z'),('PI3Ki+MEKi','PI103+PD')]
genes=['HMGCS1','HMGCR','MVD','MVK','FDPS','IDI1','SCD','SREBF1','SQLE','LDLR','INSIG1','ASNS','SREBF2','SCAP']
rows=[]
res={}
for n,f in conds:
    df=pd.read_csv(base+f'out/de_{f}.tsv',sep='\t',index_col=0)
    res[n]=df
print(res['PI3Ki'].columns.tolist())
out=[]
for g in genes:
    r=[g]
    for n,_ in conds:
        df=res[n]
        if g in df.index:
            x=df.loc[g]
            r.append(f"{x.log2FoldChange:+.2f}/{x.lfc_shrunk:+.2f}{'*' if x.padj<0.05 else ''}")
        else: r.append('NA')
    out.append(r)
print(pd.DataFrame(out,columns=['gene']+[c for c,_ in conds]).to_string(index=False))
print('format: unshrunk/shrunk; * padj<.05')
