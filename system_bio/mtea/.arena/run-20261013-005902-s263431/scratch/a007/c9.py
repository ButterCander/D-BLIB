import pandas as pd, numpy as np
P = r"D:/Cander/School/3-2/system_bio/mtea/"
x4 = pd.ExcelFile(P + 'si/MOESM4.xlsx')
rows = []
for s in x4.sheet_names:
    d = pd.read_excel(x4, s, header=1).iloc[:, :5]; d.columns = ['id', 'desc', 'score', 'p', 'sig']
    d['tid'] = 'X' + d.id.astype(int).astype(str); d = d.set_index('tid')
    for t in ('X138', 'X140', 'X141', 'X142'):
        rows.append((s, t, d.loc[t, 'desc'][:42], round(d.loc[t, 'score'], 2), d.loc[t, 'p'], bool(d.loc[t, 'sig'])))
print(pd.DataFrame(rows, columns=['sheet', 'task', 'desc', 'score', 'p', 'sig']).to_string(index=False))
# which tasks are significant (p<0.025) in all 3 isoprenoid per MEKi sheet
