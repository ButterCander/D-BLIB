import pandas as pd, numpy as np, pickle
from scipy.stats import false_discovery_control as bh, pearsonr, spearmanr
m4 = pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a008/m4.pkl','rb'))
cond = {'TAKi':('5Z','TAKi'),'MEKi':('PD0325901','MEKi'),'PI3Ki':('PI103','PI3Ki'),'PI3Ki+TAKi':('PI103+5Z','PI3Ki-TAKi'),'PI3Ki+MEKi':('PI103+PD','PI3Ki-MEKi')}
for c,(f,sh) in cond.items():
    a = m4[sh].copy(); a['task_id']='X'+a.id.astype(int).astype(str)
    for tag in ('tide_all_p10000','tide_p10000','tidee_all_p10000'):
        d = pd.read_csv(f'out/{tag}_{f}.tsv', sep='\t')
        n025 = (d.pvalue<0.025).sum(); q = bh(np.clip(2*d.pvalue.values,1e-12,1)); nb=(q<.05).sum()
        j = a.merge(d, on='task_id')
        r = pearsonr(j.score_x, j.score_y)[0] if tag.startswith('tide_') else np.nan
        print(f'{c:11} {tag:18} p<.025 {n025:3d}  BH<.05 {nb:3d}  n={len(d)}  r(authors,mine)={r:.3f}' )
