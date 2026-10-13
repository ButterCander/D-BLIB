import pandas as pd, numpy as np
P = r"D:/Cander/School/3-2/system_bio/mtea/"
L = ['CAPAN2', 'CFPAC1', 'MiaPaca2', 'PANC1']
SREBP = ['HMGCS1', 'HMGCR', 'MVK', 'PMVK', 'MVD', 'IDI1', 'FDPS', 'FDFT1', 'SQLE', 'LSS', 'CYP51A1', 'MSMO1', 'NSDHL', 'HSD17B7', 'DHCR7', 'DHCR24', 'SC5D', 'EBP', 'TM7SF2', 'ACAT2', 'INSIG1', 'LDLR', 'SCD', 'FASN']
D = {l: pd.read_csv(P + f'out/pdac_de_{l}_h24.tsv', sep='\t').set_index('gene') for l in L}
rows = []
for g in ['HMGCS1', 'HMGCR', 'MVK', 'PMVK', 'MVD', 'IDI1', 'FDPS', 'SQLE', 'LSS', 'SCD', 'INSIG1', 'LDLR', 'SREBF2', 'SCAP', 'FASN']:
    r = [g]
    for l in L:
        d = D[l]
        if g in d.index:
            star = '*' if d.loc[g, 'adj.P.Val'] < 0.05 else ''
            r.append("%+.2f%s" % (d.loc[g, 'logFC'], star))
        else:
            r.append('NA')
    rows.append(r)
print(pd.DataFrame(rows, columns=['gene'] + L).to_string(index=False))
print('program mean (24 listed genes)')
for l in L:
    d = D[l]['logFC']
    sub = d.reindex([s for s in SREBP if s in d.index])
    print(l, round(sub.mean(), 3), 'n present', len(sub), 'n up(>0)', int((sub > 0).sum()))
for l in L:
    d = D[l]['logFC']
    print(l, 'trunk4 HMGCS1/HMGCR/MVK/MVD mean', round(d.reindex(['HMGCS1', 'HMGCR', 'MVK', 'MVD']).mean(), 3))
print('n DEGs padj<.05:', {l: int((D[l]['adj.P.Val'] < 0.05).sum()) for l in L})
