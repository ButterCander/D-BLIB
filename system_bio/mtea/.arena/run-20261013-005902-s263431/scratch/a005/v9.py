import pandas as pd, numpy as np
from scipy import stats
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
gi=pd.read_csv('tsirvouli/gene_info.tsv',sep='\t',index_col=0)
sym=pd.read_csv('pdac/statin_symbols.tsv',sep='\t',index_col=0).iloc[:,0]
d=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t').set_index('Ensembl_Gene_ID')
d.index=[i.split('.')[0] for i in d.index]
d['sym']=gi.gene_name.reindex(d.index).values
d['sym2']=sym.reindex(d.index).values
d['sym']=d.sym.fillna(d.sym2)
d=d.dropna(subset=['sym']).drop_duplicates('sym').set_index('sym').drop(columns=['sym2'])
cpm=d/d.sum()*1e6
L=np.log2(cpm+1)
for line,p in [('KP4','KP4'),('MiaPaCa-2','Mia')]:
    c=L[[f'{p}-DMSO-{i}' for i in (1,2,3)]]; t=L[[f'{p}-T100-{i}' for i in (1,2,3)]]
    lf=(t.mean(axis=1)-c.mean(axis=1))
    pm=lf.reindex([g for g in SREBP if g in lf.index]).mean()
    tp=stats.ttest_ind(t.values,c.values,axis=1).pvalue
    pp=pd.Series(tp,index=L.index)
    print(line,'program mean %.2f'%pm,'n genes',len([g for g in SREBP if g in lf.index]),{g:(round(lf[g],2),round(pp[g],3)) for g in ['HMGCS1','HMGCR','MVK','MVD','SCD','LDLR']})
    # random sets
    rng=np.random.default_rng(0); expr=lf[(cpm[[f'{p}-DMSO-{i}' for i in (1,2,3)]].mean(axis=1)>5)]
    rs=np.array([expr.sample(24,random_state=int(rng.integers(1e9))).mean() for _ in range(2000)])
    print('   percentile of program mean among random 24-gene sets: %.3f'%np.mean(rs<pm))
# RMC-6236
for line in ['PSN1','SW1990']:
    x=pd.read_csv(f'pdac/GSE315130_{line}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    L=np.log2(x[[c for c in x.columns if c.startswith(line)]]+1)
    lf=L[[f'{line}_6236_{i}' for i in (1,2,3)]].mean(axis=1)-L[[f'{line}_Con_{i}' for i in (1,2,3)]].mean(axis=1)
    gs=[g for g in SREBP if g in lf.index]
    print(line,'RMC-6236 program mean %.2f'%lf[gs].mean(),'n',len(gs),'n up',int((lf[gs]>0).sum()),{g:round(lf[g],2) for g in ['HMGCS1','HMGCR','MVD','MVK','SCD']})
    lt=L[[f'{line}_Taze_{i}' for i in (1,2,3)]].mean(axis=1)-L[[f'{line}_Con_{i}' for i in (1,2,3)]].mean(axis=1)
    print('    tazemetostat program mean %.2f'%lt[gs].mean())
# MRTX1133
x=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
print(x.columns.tolist()[:40])
