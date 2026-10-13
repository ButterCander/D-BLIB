import pandas as pd, numpy as np
P = r"D:/Cander/School/3-2/system_bio/mtea/"
files = [('PI3Ki', 'de_PI103'), ('MEKi', 'de_PD0325901'), ('TAKi', 'de_5Z'), ('PI3Ki+TAKi', 'de_PI103+5Z'), ('PI3Ki+MEKi', 'de_PI103+PD')]
D = {}
for n, f in files:
    d = pd.read_csv(P + f'out/{f}.tsv', sep='\t', index_col=0)
    D[n] = d
print(D['PI3Ki'].columns.tolist(), D['PI3Ki'].shape)
genes = ['HMGCS1', 'HMGCR', 'MVD', 'MVK', 'PMVK', 'IDI1', 'FDPS', 'SCD', 'SREBF1', 'SREBF2', 'INSIG1', 'LDLR', 'SQLE', 'LSS', 'ASNS']
rows = []
for g in genes:
    r = [g]
    for n, _ in files:
        d = D[n]
        if g in d.index:
            star = '*' if d.loc[g, 'padj'] < 0.05 else ''
            r.append("%+.2f%s" % (d.loc[g, 'log2FoldChange'], star))
        else:
            r.append('NA')
    rows.append(r)
print(pd.DataFrame(rows, columns=['gene'] + [n for n, _ in files]).to_string(index=False))
# authors' MOESM2 DEG lists (apeglm LFC) for the same genes
x = pd.ExcelFile(P + 'si/MOESM2.xlsx')
for s in x.sheet_names:
    d = pd.read_excel(x, s)
    d = d.set_index('geneSymbol')
    out = {g: (round(d.loc[g, 'log2FoldChange'], 2) if g in d.index else None) for g in ['HMGCS1', 'HMGCR', 'MVD', 'MVK', 'SCD', 'SREBF1']}
    print('MOESM2', s, len(d), out)
# Tsirvouli official
off = pd.read_csv(P + 'tsirvouli_official/DEGs_perContrast_all.csv')
gi = pd.read_csv(P + 'tsirvouli/gene_info.tsv', sep='\t', index_col=0)
off['sym'] = gi.reindex(off.Gene_name)['gene_name'].values
print(off.columns.tolist())
print(sorted(off.Contrast.unique())[:60])
for c in ['PI', 'PD', '5Z', 'PIPD', 'PI5Z']:
    sub = off[off.Contrast == f'Diff_{c}_24h']
    r = {}
    for g in ['HMGCS1', 'HMGCR', 'MVK', 'MVD', 'IDI1', 'SCD', 'INSIG1', 'LDLR']:
        t = sub[sub.sym == g]
        r[g] = (round(t.logFC.iloc[0], 2), float('%.3g' % t.adj_P_Val.iloc[0]) if 'adj_P_Val' in t.columns else None) if len(t) else None
    print(c, r)
