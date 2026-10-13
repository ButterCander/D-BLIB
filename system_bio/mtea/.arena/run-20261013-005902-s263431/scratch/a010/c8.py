import pandas as pd, numpy as np, warnings
from scipy import stats
warnings.filterwarnings('ignore')
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
def report(name,lf,p=None):
    s=lf.reindex([g for g in SREBP if g in lf.index]).dropna()
    t=lf.reindex([g for g in TR if g in lf.index])
    rng=np.random.default_rng(1); pool=lf.dropna().values
    rnd=np.array([rng.choice(pool,len(s),replace=False).mean() for _ in range(3000)])
    print(f"{name:34} n={len(s):2d} up={int((s>0).sum()):2d} prog={s.mean():+.2f} frac_rand>=obs={(rnd>=s.mean()).mean():.4f} | "+" ".join(f"{g}:{t[g]:+.2f}" for g in TR if g in t.index and pd.notna(t[g])))
# GSE206267
c=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
print(list(c.columns))
gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)['gene_name']
c=c[c.sum(axis=1)>=30]; c['sym']=gi.reindex(c.index).values; c=c.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym')
nc=c/c.sum()*c.sum().median(); L=np.log2(nc+1)
for line,pre in (('KP4','KP4'),('MiaPaCa2','Mia')):
    a=L[[f'{pre}-T100-{i}' for i in (1,2,3)]]; b=L[[f'{pre}-DMSO-{i}' for i in (1,2,3)]]
    lf=a.mean(axis=1)-b.mean(axis=1); report('GSE206267 '+line,lf)
# GSE315130
for line in ('PSN1','SW1990'):
    d=pd.read_csv(f'pdac/GSE315130_{line}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t')
    print(line,[c for c in d.columns][:20])
    d=d.drop_duplicates('Symbol').set_index('Symbol')
    for tag in ('6236','Taze','Combo'):
        cols=[f'{line}_{tag}_{i}' for i in (1,2,3)]
        if not all(c in d.columns for c in cols): continue
        lf=np.log2(d[cols]+1).mean(axis=1)-np.log2(d[[f'{line}_Con_{i}' for i in (1,2,3)]]+1).mean(axis=1)
        report(f'GSE315130 {line} {tag}',lf)
# GSE279749
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t')
print(d.columns.tolist()[:40], d.shape)
