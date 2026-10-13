import numpy as np, pandas as pd
from scipy import stats
cnt = pd.read_csv('GSE285616_counts.csv.gz', index_col='Symbol')
cnt = cnt[~cnt.index.duplicated()]
print(cnt.shape)
# median-of-ratios size factors
x = cnt.loc[(cnt>0).all(axis=1)]
lg = np.log(x)
sf = np.exp((lg.sub(lg.mean(axis=1), axis=0)).median(axis=0))
print(sf.round(2).tolist())
N = cnt/sf
conds = ['control','PI103','PD0325901','5Z','PI103+5Z','PI103+PD']
cols = list(cnt.columns)
idx = {c:[cols[b*6+i] for b in range(4)] for i,c in enumerate(conds)}
L = np.log2(N+1)  # log2(norm+1) for paired t
res = {}
keep = cnt.sum(axis=1)>=10
for c in conds[1:]:
    t = N[idx[c]].values; k = N[idx['control']].values
    # paired within batch
    lfc_ratio = np.log2((t.mean(axis=1)+0.5)/(k.mean(axis=1)+0.5))
    d = np.log2(t+1)-np.log2(k+1)
    tt = stats.ttest_1samp(d, 0, axis=1)
    p = pd.Series(tt.pvalue, index=cnt.index)
    from scipy.stats import false_discovery_control as bh
    pp = p.copy(); ok = keep & p.notna()
    q = pd.Series(np.nan, index=cnt.index); q[ok] = bh(p[ok].values)
    res[c] = pd.DataFrame({'lfc':lfc_ratio, 'lfc_paired':d.mean(axis=1), 'q':q, 'mean_ctrl':k.mean(axis=1)}, index=cnt.index)
import pickle; pickle.dump(res, open('.arena/run-20261013-005902-s263431/scratch/a004/de_simple.pkl','wb'))
genes = ['HMGCS1','HMGCR','MVD','MVK','SCD','SREBF1','IDI1','FDPS','INSIG1','LDLR','SQLE','ASNS']
for c in conds[1:]:
    print(c); print(res[c].loc[genes].round(3).to_string())
