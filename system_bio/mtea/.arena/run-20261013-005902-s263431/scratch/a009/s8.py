import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
files={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki+TAKi':'PI103+5Z','PI3Ki+MEKi':'PI103+PD'}
D={k:pd.read_csv(f'out/de_{v}.tsv',sep='\t') for k,v in files.items()}
print(D['MEKi'].columns.tolist(), len(D['MEKi']))
c0=D['MEKi'].columns
idcol=c0[0]
for k in D: D[k]=D[k].set_index(D[k].columns[0])
print(D['MEKi'].head(2))
for gn in ['HMGCS1','HMGCR','MVK','MVD','SCD','SREBF1','IDI1','FDPS','ASNS']:
    print(f'{gn:7s}', ' '.join(f"{k}:{D[k].loc[gn,'lfc_shrunk']:+.2f}/{D[k].loc[gn,'log2FoldChange']:+.2f}({D[k].loc[gn,'padj']:.3f})" for k in D))
# DEG counts |lfc|>=1 & padj<.05
for k in D:
    d=D[k]; print(k,'DEG(|lfc_shrunk|>=1)',int(((d.padj<0.05)&(d.lfc_shrunk.abs()>=1)).sum()),'(|lfc_unshrunk|>=1)',int(((d.padj<0.05)&(d.log2FoldChange.abs()>=1)).sum()),'up',int(((d.padj<0.05)&(d.lfc_shrunk>=1)).sum()),'down',int(((d.padj<0.05)&(d.lfc_shrunk<=-1)).sum()))
