import pandas as pd, numpy as np
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TRUNK=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
def prog(d,tr,ct,genes=SREBP):
    L=np.log2(d[tr+ct].astype(float)+1)
    lfc=L[tr].mean(axis=1)-L[ct].mean(axis=1)
    x=lfc.reindex(genes).dropna()
    return x.mean(), len(x), int((x>0).sum()), lfc
for line in ['PSN1','SW1990']:
    d=pd.read_csv(f'pdac/GSE315130_{line}_cell_line_gene_sample_FPKM_with_symbol.txt.gz',sep='\t').drop_duplicates('Symbol').set_index('Symbol')
    for lab,k in [('RMC-6236','6236'),('Taze','Taze'),('Combo','Combo')]:
        m,n,up,lfc=prog(d,[f'{line}_{k}_{i}' for i in (1,2,3)],[f'{line}_Con_{i}' for i in (1,2,3)])
        print(line,lab,'program mean %.2f (n=%d, up %d)'%(m,n,up),'trunk %.2f'%lfc.reindex(TRUNK).mean(), ' HMGCS1 %.2f HMGCR %.2f'%(lfc['HMGCS1'],lfc['HMGCR']))
d=pd.read_csv('pdac/GSE279749_gene_fpkm_cell_lines.txt.gz',sep='\t').drop_duplicates('gene_name').set_index('gene_name')
for lab,k in [('AsPC-1 MRTX1133','A'),('PANC-1 MRTX1133','P'),('CFPAC-1 BI-2865','C'),('SW1990 MRTX1133','S9')]:
    m,n,up,lfc=prog(d,[f'{k}_M_{i}' for i in (1,2,3)],[f'{k}_NT_{i}' for i in (1,2,3)])
    print(lab,'program mean %.2f (n=%d, up %d)'%(m,n,up),'trunk %.2f'%lfc.reindex(TRUNK).mean(), 'HMGCS1 %.2f HMGCR %.2f MVD %.2f IDI1 %.2f'%(lfc['HMGCS1'],lfc['HMGCR'],lfc['MVD'],lfc['IDI1']))
