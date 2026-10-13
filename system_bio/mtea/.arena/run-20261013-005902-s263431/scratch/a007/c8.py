import pandas as pd, numpy as np, json
from scipy.stats import false_discovery_control as bh
P = r"D:/Cander/School/3-2/system_bio/mtea/"
x4 = pd.ExcelFile(P + 'si/MOESM4.xlsx'); x5 = pd.ExcelFile(P + 'si/MOESM5.xlsx')
S = {}
for s in x4.sheet_names:
    d = pd.read_excel(x4, s, header=1).iloc[:, :5]; d.columns = ['id', 'desc', 'score', 'p', 'sig']
    d['tid'] = 'X' + d.id.astype(int).astype(str); S[s] = d.set_index('tid')
for comb, singles, a, b in (('PI3Ki-TAKi', ('PI3Ki', 'TAKi'), 'PI3Ki', 'TAKi'), ('PI3Ki-MEKi', ('PI3Ki', 'MEKi'), 'PI3Ki', 'MEKi')):
    m = pd.read_excel(x5, comb)
    m.columns = ['task_id', 'description', 'A', 'B', 'AB', 'direction'] + list(m.columns[6:])
    q = pd.Series(bh(np.clip(2 * S[comb].p.values, 1e-12, 1)), index=S[comb].index)
    rows = []
    for r in m.itertuples():
        t = r.task_id
        pa, pb, pab = S[a].loc[t, 'p'], S[b].loc[t, 'p'], S[comb].loc[t, 'p']
        sa, sb, sab = S[a].loc[t, 'score'], S[b].loc[t, 'score'], S[comb].loc[t, 'score']
        strong = (r.A == 0) and (r.B == 0)
        rows.append((t, r.description[:40], r.direction, 'strong' if strong else 'weak', round(sa, 3), round(pa, 3), round(sb, 3), round(pb, 3), round(sab, 3), round(pab, 4), round(q[t], 3),
                     abs(sab) - max(abs(sa), abs(sb))))
    df = pd.DataFrame(rows, columns=['t', 'desc', 'dir', 'type', 'sA', 'pA', 'sB', 'pB', 'sAB', 'pAB', 'qAB(2p)', 'margin_vs_best_single_abs'])
    print(comb, len(df))
    st = df[df.type == 'strong']
    print(st.drop(columns=['margin_vs_best_single_abs']).to_string(index=False))
    wk = df[df.type == 'weak']
    print('weak n', len(wk), 'median margin (|AB| - max |single| using sheet scores) %.3f' % wk.margin_vs_best_single_abs.median(), 'frac <0.10: %.2f' % (wk.margin_vs_best_single_abs < 0.10).mean(),
          'p of AB<0.025 count', int((wk.pAB < 0.025).sum()), 'q<.05 count', int((wk['qAB(2p)'] < .05).sum()))
    # for weak: are the single drugs' scores significant in the sheet?
    print('weak: share where a single-drug p<0.025 :', round(((wk.pA < 0.025) | (wk.pB < 0.025)).mean(), 2))
    print('overall: calls with AB p<.025: ', int((df.pAB < .025).sum()), ' q<.05:', int((df['qAB(2p)'] < .05).sum()), ' q<.10:', int((df['qAB(2p)'] < .10).sum()))
# slide-script numbers
for s in ['PI3Ki', 'PI3Ki-TAKi', 'PI3Ki-MEKi']:
    d = S[s]; sig = d[d.p < 0.025]
    print(s, 'sig', len(sig), 'down', int((sig.score < 0).sum()), 'frac down %.2f' % (sig.score < 0).mean(), 'median |score| of sig %.2f' % sig.score.abs().median())
tot = pd.concat([S[s][S[s].p < 0.025] for s in S]); print('all five: sig', len(tot), 'down frac %.3f' % (tot.score < 0).mean())
# asn
print({s: round(S[s].loc['X76', 'score'], 2) for s in S}, {s: S[s].loc['X76', 'p'] for s in S})
tg = json.load(open(P + 'out/task_genes.json')); print('X76 genes', tg.get('X76'))
