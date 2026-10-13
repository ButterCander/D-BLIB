import json, numpy as np, pandas as pd
from scipy import stats
P = r"D:/Cander/School/3-2/system_bio/mtea/"
def load(p):
    return pd.read_csv(P + p, sep="\t").set_index("task_id")["score"]
GS = {"PI": "PI103", "PD": "PD0325901", "OXO": "5Z", "PIOXO": "PI103+5Z", "PIPD": "PI103+PD"}
for lab, gp, tp in (("TIDE", "out/tide_p10000_%s.tsv", "out/tidex_tsir_%s_24h.tsv"), ("TIDE-ess", "out/tidee_%s.tsv", "out/tideex_tsir_%s_24h.tsv")):
    G = {k: load(gp % v) for k, v in GS.items()}
    T = {k: load(tp % k) for k in GS}
    M = np.array([[stats.spearmanr(*pd.concat([G[r], T[c]], axis=1).dropna().values.T)[0] for c in GS] for r in GS])
    d = np.diag(M); off = M[~np.eye(5, dtype=bool)]
    print(lab, 'matched', np.round(d, 2), 'mean %.2f' % d.mean(), '| unmatched mean %.2f max %.2f min %.2f' % (off.mean(), off.max(), off.min()))
# naive mean vs TIDE
tasks = json.load(open(P + "out/task_genes.json"))
for c in ["PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]:
    de = pd.read_csv(P + f"out/de_{c}.tsv", sep="\t", index_col=0)
    lf = de["lfc_shrunk"].where(de["padj"].fillna(1) < .05, 0.0)
    naive = pd.Series({t: lf.reindex(g).dropna().mean() for t, g in tasks.items() if lf.reindex(g).dropna().size >= 1})
    tide = pd.read_csv(P + f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")["score"]
    j = pd.concat([tide, naive], axis=1, keys=["tide", "naive"]).dropna()
    print(c, 'n', len(j), 'pearson %.2f spearman %.2f' % (stats.pearsonr(j.tide, j.naive)[0], stats.spearmanr(j.tide, j.naive)[0]))
# TIDE vs authors scores
x = pd.ExcelFile(P + 'si/MOESM4.xlsx')
m = {'PI3Ki': 'PI103', 'MEKi': 'PD0325901', 'TAKi': '5Z', 'PI3Ki-TAKi': 'PI103+5Z', 'PI3Ki-MEKi': 'PI103+PD'}
for s, c in m.items():
    d = pd.read_excel(x, s, header=1).iloc[:, :5]; d.columns = ['id', 'desc', 'score', 'p', 'sig']
    d['task_id'] = 'X' + d.id.astype(int).astype(str); d = d.set_index('task_id')
    mine = pd.read_csv(P + f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")
    j = pd.concat([d.score, mine.score], axis=1, keys=['auth', 'mine']).dropna()
    print(s, '(sheet name) vs', c, 'n', len(j), 'pearson %.2f spearman %.2f' % (stats.pearsonr(j.auth, j.mine)[0], stats.spearmanr(j.auth, j.mine)[0]))
# cross: sheet TAKi vs my PD (to test the header swap)
for s, c in (('TAKi', 'PD0325901'), ('MEKi', '5Z')):
    d = pd.read_excel(x, s, header=1).iloc[:, :5]; d.columns = ['id', 'desc', 'score', 'p', 'sig']
    d['task_id'] = 'X' + d.id.astype(int).astype(str); d = d.set_index('task_id')
    mine = pd.read_csv(P + f"out/tide_p10000_{c}.tsv", sep="\t").set_index("task_id")
    j = pd.concat([d.score, mine.score], axis=1, keys=['auth', 'mine']).dropna()
    print('swapped test', s, 'vs', c, 'pearson %.2f' % stats.pearsonr(j.auth, j.mine)[0])
