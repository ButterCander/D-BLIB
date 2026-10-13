import pandas as pd, numpy as np
gs=pd.read_csv('external/crispr_gene_summary.tsv',sep='\t'); gs['sym']=gs.iloc[:,0].str.split(' ').str[0]; gs=gs.set_index('sym')
for g in ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','GGPS1','FNTA','SQLE','LSS','SCD','SCAP','RRM1','RRM2','KRAS','ASNS','PHGDH','ODC1']:
    if g in gs.index: r=gs.loc[g]; print(f'{g:7} all-lines mean {r["mean"]:.2f} frac_dep {r.frac_dep:.2f} n {int(r.n)}')
print('genes with frac_dep>=0.95:',(gs.frac_dep>=0.95).sum(),'of',len(gs))
gr=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0)
print(gr.shape, [c for c in gr.columns if c.startswith('HMGCS1')][:3], list(gr.columns[:3]))
