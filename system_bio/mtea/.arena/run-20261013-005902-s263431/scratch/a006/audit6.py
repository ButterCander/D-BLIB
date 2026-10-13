import pandas as pd, numpy as np, json, glob
B='D:/Cander/School/3-2/system_bio/mtea/'
print('statin positive control: rank (by score desc, 187 tasks) / score / p for X138 X140 X141')
for line in ['PANC1','MiaPaca2']:
    for dose in ['High','Low']:
        for st in ['Fluvastatin','Lovastatin','Simvastatin']:
            d=pd.read_csv(B+'out/tidex_statin_%s_%s_%s.tsv'%(line,st,dose),sep='\t').set_index('task_id')
            d['rk']=d.score.rank(ascending=False,method='min')
            print(line,dose,st,' '.join('%s rk%d sc%.2f p%.3f'%(t,d.rk[t],d.score[t],d.pvalue[t]) for t in ['X138','X140','X141']),'n_p<.025=%d'%(d.pvalue<.025).sum())
g=json.load(open(B+'out/eflux_grid.json'))
vals=[v for x in g for v in x['pred'].values()]
print('E-Flux grid settings',len(g),'n values',len(vals),'min %.1f max %.1f'%(min(vals),max(vals)))
print('settings with rho:',[round(x['rho'],2) for x in g][:20])
mins=[min(x['pred'].values()) for x in g]; maxs=[max(x['pred'].values()) for x in g]
print('per-setting min range',min(mins),max(mins),'max range',min(maxs),max(maxs))
