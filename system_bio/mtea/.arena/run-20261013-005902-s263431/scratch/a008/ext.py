import pandas as pd, numpy as np
SREBP = ['HMGCS1','HMGCR','MVK','PMVK','MVD','IDI1','FDPS','FDFT1','SQLE','LSS','CYP51A1','MSMO1','NSDHL','HSD17B7','DHCR7','DHCR24','SC5D','EBP','TM7SF2','ACAT2','INSIG1','LDLR','SCD','FASN']
TR=['HMGCS1','HMGCR','MVK','MVD']
def rep(name, lf):
    p=lf.reindex([s for s in SREBP if s in lf.index])
    print(f'{name:28} program mean {p.mean():+.2f} (n={p.notna().sum()}, up {int((p>0).sum())})  ' + ' '.join(f'{g} {lf.get(g,np.nan):+.2f}' for g in TR))
# GSE279749 FPKM
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
for pre,nm in [('A','AsPC-1'),('P','PANC-1'),('C','CFPAC-1'),('S9','SW1990')]:
    nt=np.log2(d[[f'{pre}_NT_{i}' for i in (1,2,3)]]+1).mean(axis=1); m=np.log2(d[[f'{pre}_M_{i}' for i in (1,2,3)]]+1).mean(axis=1)
    rep('GSE279749 '+nm, m-nt)
# GSE315130
for f in ['PSN1','SW1990']:
    h=pd.read_csv(f'pdac/GSE315130_{f}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    con=np.log2(h[[f'{f}_Con_{i}' for i in (1,2,3)]]+1).mean(axis=1)
    for k in ('6236','Combo','Taze'):
        t=np.log2(h[[f'{f}_{k}_{i}' for i in (1,2,3)]]+1).mean(axis=1)
        rep(f'GSE315130 {f} {k}', t-con)
# GSE206267 raw counts: simple CPM log2
c=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
import re
# map ensembl -> symbol via Human-GEM genes or statin_symbols
try:
    sm=pd.read_csv('pdac/statin_symbols.tsv',sep='\t'); print(sm.head(2))
except Exception as e: print(e)
sm=sm.drop_duplicates('gene_id').set_index('gene_id').gene_short_name
c=c.loc[c.index.isin(sm.index)]
c.index=sm.reindex(c.index).values
c=c.groupby(level=0).sum()
cpm=np.log2(c/c.sum()*1e6+1)
for pre,nm in [('KP4','KP4'),('Mia','MiaPaCa-2')]:
    dm=cpm[[f'{pre}-DMSO-{i}' for i in (1,2,3)]]; tr=cpm[[f'{pre}-T100-{i}' for i in (1,2,3)]]
    lf=tr.mean(axis=1)-dm.mean(axis=1)
    rep('GSE206267 '+nm+' 48h (logCPM)', lf)
