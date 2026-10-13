import pandas as pd, numpy as np, json, gzip, openpyxl
from scipy import stats
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
print("=== E. PDAC program mean (24 listed genes) from normalised counts, log2(x+1) mean diff vs NT ===")
for line in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    d=pd.read_csv(f'pdac/GSE238198_Pdac.{line}.genes.normalizedCounts.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    L=np.log2(d+1)
    nt=[c for c in d.columns if '_NT_' in c]
    out={}
    for h in ['3h','12h','24h']:
        tr=[c for c in d.columns if f'_{h}_' in c]
        lf=L[tr].mean(axis=1)-L[nt].mean(axis=1)
        g=[s for s in SREBP if s in lf.index]
        out[h]=round(lf[g].mean(),2)
        if h=='24h':
            out['n_genes']=len(g)
            out['HMGCR24']=round(lf['HMGCR'],2); out['HMGCS1_24']=round(lf['HMGCS1'],2); out['PMVK24']=round(lf['PMVK'],2)
            out['frac_up']=round((lf[g]>0).mean(),2)
    print(line,out)
print()
print("=== GSE206267 KP4 / MiaPaCa-2 trametinib 100 nM 48h (raw counts, my CPM-log2 diff) ===")
d=pd.read_csv('pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz',sep='\t',index_col=0)
print(d.shape)
sym=pd.read_csv('mteapy/src/mteapy/data/HumanGEM_genes.tsv',sep='\t')
m=dict(zip(sym.geneID,sym.geneSymbol))
d=d[d.index.isin(m)]; d.index=[m[i] for i in d.index]; d=d[~d.index.duplicated()]
cpm=d/d.sum()*1e6; L=np.log2(cpm+1)
for line,a,b in [('KP4',['KP4-DMSO-1','KP4-DMSO-2','KP4-DMSO-3'],['KP4-T100-1','KP4-T100-2','KP4-T100-3']),('Mia',['Mia-DMSO-1','Mia-DMSO-2','Mia-DMSO-3'],['Mia-T100-1','Mia-T100-2','Mia-T100-3'])]:
    lf=L[b].mean(axis=1)-L[a].mean(axis=1)
    g=[s for s in SREBP if s in lf.index]
    print(line,{k:round(lf[k],2) for k in ['HMGCS1','HMGCR','MVK','MVD']}, 'program',round(lf[g].mean(),2),len(g))
print()
print("=== GSE315130 PSN1/SW1990 RMC-6236 vs Con (FPKM log2(x+1)) ===")
for nm in ['PSN1','SW1990']:
    d=pd.read_csv(f'pdac/GSE315130_{nm}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    cols=d.columns
    a=[c for c in cols if '_Con_' in c]; b=[c for c in cols if '_6236_' in c]
    L=np.log2(d[a+b].astype(float)+1)
    lf=L[b].mean(axis=1)-L[a].mean(axis=1)
    g=[s for s in SREBP if s in lf.index]
    print(nm,a,b,{k:round(lf[k],2) for k in ['HMGCS1','HMGCR','MVK','MVD']},'program',round(lf[g].mean(),2),len(g))
print()
print("=== GSE279749 KRAS inhibitors 24 h ===")
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t')
print(d.columns.tolist()[:40])
