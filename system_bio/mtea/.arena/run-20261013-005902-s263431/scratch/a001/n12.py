import pandas as pd, numpy as np
print('line statin dose: ranks (of 187 by score desc) X138 chol, X140 mevalonate, X141 FPP ; p')
for line in ['PANC1','MiaPaca2']:
    for st in ['Fluvastatin','Lovastatin','Simvastatin']:
        for lev in ['High','Low']:
            d=pd.read_csv(f'out/tidex_statin_{line}_{st}_{lev}.tsv',sep='\t').set_index('task_id')
            r=d.score.rank(ascending=False,method='min')
            s=[(int(r[t]),round(d.score[t],2),d.pvalue[t]) for t in ('X138','X140','X141')]
            print(line,st,lev,s, 'n=',d.score.notna().sum())
