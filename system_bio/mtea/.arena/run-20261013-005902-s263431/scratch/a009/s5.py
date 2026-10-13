import pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
p=pd.read_csv('external/crispr_pancreas_rows.tsv',sep='\t',index_col=0)
g=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0)
s=pd.read_csv('external/crispr_gene_summary.tsv',sep='\t',index_col=0)
for d in (p,g): d.columns=[c.split(' (')[0] for c in d.columns]
s.index=[i.split(' (')[0] for i in s.index]
print('pancreas rows',p.shape,'gastric rows',g.shape, 'summary n genes',len(s))
def row(gn):
    v=p[gn].dropna(); vg=g[gn].dropna()
    return f'{gn:7s} PDAC n={len(v)} median={v.median():.2f} frac<-0.5={np.mean(v<-0.5):.2f} | gastric n={len(vg)} median={vg.median():.2f} frac={np.mean(vg<-0.5):.2f} | allCL frac_dep={s.loc[gn,"frac_dep"]:.2f} mean={s.loc[gn,"mean"]:.2f} n={s.loc[gn,"n"]:.0f}'
for gn in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','FNTB','SQLE','LSS','DHCR7','FDFT1','SCD','SCAP','MBTPS1','SREBF2','RRM1','RRM2','CTPS1','ASNS','LDLR','INSIG1']:
    if gn in p.columns: print(row(gn))
    else: print(gn,'missing')
print(p.index.tolist())
