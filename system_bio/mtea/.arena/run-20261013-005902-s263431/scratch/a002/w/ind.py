import pandas as pd, numpy as np
from scipy import stats
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
def summ(l,name,pc=0.5):
    # l: DataFrame gene x (ctrl cols, trt cols) -> lfc series
    return l
res={}
# GSE238198 24h
for ln in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    d=pd.read_csv(f'pdac/GSE238198_Pdac.{ln}.genes.normalizedCounts.txt.gz',sep='\t',index_col=0)
    nt=[c for c in d.columns if '_NT_' in c]; t=[c for c in d.columns if '_24h_' in c]
    lfc=np.log2((d[t].mean(axis=1)+0.5)/(d[nt].mean(axis=1)+0.5))
    res[('tram24',ln)]=lfc
# GSE206267 48h raw counts -> need symbol; map via gene id? check if symbol available elsewhere
d=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
# map ensembl->symbol using Human-GEM gene file or DepMap? try mteapy data HumanGEM_genes.tsv (geneID, geneSymbol)
hg=pd.read_csv('mteapy/src/mteapy/data/HumanGEM_genes.tsv',sep='\t').drop_duplicates('geneID').set_index('geneID')['geneSymbol']
lib=d.sum()
cpm=d/lib*1e6
for cl,pre in (('KP4','KP4'),('MiaPaca2','Mia')):
    ctrl=[c for c in d.columns if c.startswith(pre+'-DMSO')]; tr=[c for c in d.columns if c.startswith(pre+'-T100')]
    lfc=np.log2((cpm[tr].mean(axis=1)+0.5)/(cpm[ctrl].mean(axis=1)+0.5))
    # filter low expression
    keep=cpm[ctrl+tr].mean(axis=1)>2
    lfc=lfc[keep]
    lfc.index=[hg.get(i,None) for i in lfc.index]; lfc=lfc[lfc.index.notna()]; lfc=lfc[~lfc.index.duplicated()]
    res[('tram48',cl)]=lfc
# GSE279749
e=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
for k,nm in (('A','AsPC1'),('P','PANC1'),('C','CFPAC1'),('S9','SW1990')):
    nt=[f'{k}_NT_{i}' for i in (1,2,3)]; m=[f'{k}_M_{i}' for i in (1,2,3)]
    lfc=np.log2((e[m].mean(axis=1)+0.5)/(e[nt].mean(axis=1)+0.5)); res[('kras24',nm)]=lfc
# GSE315130
for ln in ('PSN1','SW1990'):
    g=pd.read_csv(f'pdac/GSE315130_{ln}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    lfc=np.log2((g[[f'{ln}_6236_{i}' for i in (1,2,3)]].mean(axis=1)+0.5)/(g[[f'{ln}_Con_{i}' for i in (1,2,3)]].mean(axis=1)+0.5))
    res[('rmc6236',ln)]=lfc
print(f"{'dataset':18}{'prog24':>8}{'n':>4}{'up':>4}{'trunk6':>8}{'HMGCS1':>8}{'HMGCR':>8}{'MVK':>7}{'MVD':>7}")
for k,l in res.items():
    s=l.reindex([g for g in SREBP if g in l.index]).dropna()
    t=l.reindex(TR)
    print(f"{k[0]+'/'+k[1]:18}{s.mean():8.2f}{len(s):4d}{(s>0).sum():4d}{t.mean():8.2f}{l.get('HMGCS1',np.nan):8.2f}{l.get('HMGCR',np.nan):8.2f}{l.get('MVK',np.nan):7.2f}{l.get('MVD',np.nan):7.2f}")
import pickle; pickle.dump(res,open('.arena/run-20261013-005902-s263431/scratch/a002/w/ind_res.pkl','wb'))
