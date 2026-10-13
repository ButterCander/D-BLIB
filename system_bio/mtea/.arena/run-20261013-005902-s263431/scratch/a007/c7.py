import pandas as pd, numpy as np
from scipy import stats
P = r"D:/Cander/School/3-2/system_bio/mtea/"
x = pd.ExcelFile(P + 'si/MOESM2.xlsx')
mine = {c: pd.read_csv(P + f'out/de_{c}.tsv', sep='\t', index_col=0)['lfc_shrunk'] for c in ['PD0325901', '5Z', 'PI103', 'PI103+5Z', 'PI103+PD']}
print(pd.read_csv(P + 'out/de_PI103.tsv', sep='\t', index_col=0).head(2).to_string())
sym_ok = list(mine['PI103'].index[:5]); print(sym_ok)
for s in x.sheet_names:
    d = pd.read_excel(x, s).drop_duplicates('geneSymbol').set_index('geneSymbol')['log2FoldChange']
    out = {}
    for c, v in mine.items():
        j = pd.concat([d, v], axis=1, join='inner').dropna()
        out[c] = (len(j), round(stats.pearsonr(j.iloc[:, 0], j.iloc[:, 1])[0], 2))
    print(s, len(d), out)
