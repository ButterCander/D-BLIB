import sys, pandas as pd, numpy as np, itertools
sys.stdout.reconfigure(encoding='utf-8')
base='D:/Cander/School/3-2/system_bio/mtea/'
print('statin positive control: rank (1=highest score among 187 tasks), p')
for line in ['PANC1','MiaPaca2']:
    for dose in ['High','Low']:
        for st in ['Fluvastatin','Lovastatin','Simvastatin']:
            f=base+f'out/tidex_statin_{line}_{st}_{dose}.tsv'
            df=pd.read_csv(f,sep='\t')
            df['rank']=df.score.rank(ascending=False,method='min')
            s=df.set_index('task_id')
            out=[]
            for t in ['X138','X140','X141']:
                r=s.loc[t]
                out.append(f"{t}: rank {int(r['rank'])} score {r.score:.2f} p {r.pvalue if 'pvalue' in r else float('nan')}")
            print(line,dose,st,' | '.join(out), '| n tasks',len(df))
