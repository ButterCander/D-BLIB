import pandas as pd, numpy as np
from scipy import stats
pd.set_option('display.width',250)
M4='si/MOESM4.xlsx'
d={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    raw=pd.read_excel(M4,sheet_name=s,header=None)
    df=raw.iloc[2:,:5].copy(); df.columns=['id','desc','score','p','sig']; df=df.dropna(subset=['id'])
    df['id']='X'+df.id.astype(int).astype(str); df['score']=pd.to_numeric(df.score); df['p']=pd.to_numeric(df.p); d[s]=df.set_index('id')
for t in ['X138','X140','X141','X76']:
    print(t, d['PI3Ki'].loc[t,'desc'] if t in d['PI3Ki'].index else '')
    for s in d: 
        if t in d[s].index: r=d[s].loc[t]; print('   %-11s score %.2f p %.4f sig=%s'%(s,r.score,r.p,r.sig))
# mapping check sheet vs my reproduction
mine={'TAKi':'5Z','MEKi':'PD0325901','PI3Ki':'PI103','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
rep={k:pd.read_csv(f'out/tide_all_p10000_{v}.tsv',sep='\t').set_index('task_id') for k,v in mine.items()}
for s in ['TAKi','MEKi']:
    for m in ['TAKi','MEKi']:
        j=d[s].join(rep[m][['score']],rsuffix='_m',how='inner')
        print('sheet',s,'vs my',m,'n',len(j),'pearson %.2f'%np.corrcoef(j.score,j.score_m)[0,1],'spearman %.2f'%stats.spearmanr(j.score,j.score_m)[0])
for s in ['PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    j=d[s].join(rep[s][['score']],rsuffix='_m',how='inner'); print('sheet',s,'n',len(j),'pearson %.2f'%np.corrcoef(j.score,j.score_m)[0,1])
# BH of my reproduction
from scipy.stats import false_discovery_control as bh
for s in rep:
    p=rep[s].pvalue.values; q=bh(np.clip(2*p,1e-12,1)); print('my repro',s,'p<.025',(p<0.025).sum(),'BH<.05',(q<.05).sum())
# fraction down among significant
for s in ['PI3Ki','PI3Ki-MEKi','PI3Ki-TAKi']:
    x=d[s][d[s].p<0.025]; print(s,len(x),'down',(x.score<0).sum(),'%.0f%%'%(100*(x.score<0).mean()))
# X76 unique gene?
