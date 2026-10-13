from h import *
mine = {"TAKi":"5Z","MEKi":"PD0325901","PI3Ki":"PI103","PI3Ki-TAKi":"PI103+5Z","PI3Ki-MEKi":"PI103+PD"}
for sheet in auth:
    A = a_scores(sheet)
    row=[]
    for cond,suf in mine.items():
        m = my(suf)
        j = A.join(m, how='inner', rsuffix='_m')
        r = np.corrcoef(j['score'], j['score_m'])[0,1]
        rs = stats.spearmanr(j['score'], j['score_m'])[0]
        row.append(f"{cond}:{r:.2f}/{rs:.2f}")
    print(f"sheet {sheet:12s}", " ".join(row))
