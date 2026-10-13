import pandas as pd, numpy as np
from scipy import stats
P = r"D:/Cander/School/3-2/system_bio/mtea/"
SREBP = ['HMGCS1', 'HMGCR', 'MVK', 'PMVK', 'MVD', 'IDI1', 'FDPS', 'FDFT1', 'SQLE', 'LSS', 'CYP51A1', 'MSMO1', 'NSDHL', 'HSD17B7', 'DHCR7', 'DHCR24', 'SC5D', 'EBP', 'TM7SF2', 'ACAT2', 'INSIG1', 'LDLR', 'SCD', 'FASN']
# GSE279749
d = pd.read_csv(P + 'pdac/GSE279749_gene_fpkm_cell_lines.txt.gz', sep='\t')
print(d.columns.tolist())
d = d.drop_duplicates('gene_name').set_index('gene_name')
cols = d.columns.tolist()
pref = sorted({c.rsplit('_', 1)[0] for c in cols if c.endswith(('_1', '_2', '_3'))})
print(pref)
res = {}
for p in pref:
    if p.endswith('_NT') or p.endswith('_M') or True:
        pass
lines = sorted({p.rsplit('_', 1)[0] for p in pref})
for l in lines:
    nt = [f'{l}_NT_{i}' for i in (1, 2, 3)]
    for tr in [c.rsplit('_', 1)[0] for c in cols if c.startswith(l + '_') and '_NT_' not in c]:
        pass
print(lines)
for l in lines:
    treat = sorted({c.rsplit('_', 1)[0] for c in cols if c.startswith(l + '_') and '_NT' not in c})
    for t in treat:
        tcols = [c for c in cols if c.startswith(t + '_')]
        ncols = [f'{l}_NT_{i}' for i in (1, 2, 3) if f'{l}_NT_{i}' in cols]
        if len(tcols) and len(ncols):
            lf = np.log2(d[tcols] + 1).mean(axis=1) - np.log2(d[ncols] + 1).mean(axis=1)
            sub = lf.reindex([s for s in SREBP if s in lf.index])
            print(l, t, 'program mean %.2f' % sub.mean(), 'n', len(sub), 'up', int((sub > 0).sum()), '| HMGCS1 %.2f HMGCR %.2f MVK %.2f MVD %.2f' % tuple(lf.reindex(['HMGCS1', 'HMGCR', 'MVK', 'MVD'])))
