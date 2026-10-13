import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh
P='si/MOESM4.xlsx'
S={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel(P,sheet_name=s,header=None).iloc[2:,[0,1,2,3]]; d.columns=['id','desc','score','p']
    d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p); S[s]=d.dropna(subset=['p']).reset_index(drop=True)
    q=bh(np.clip(2*S[s].p.values,1e-12,1))
    print(s,len(S[s]),'p<.025:',(S[s].p<0.025).sum(),'up',((S[s].p<0.025)&(S[s].score>0)).sum(),'down',((S[s].p<0.025)&(S[s].score<0)).sum(),'BH<.05:',(q<.05).sum())
# sheet identity: correlate with my TIDE of 5Z and PD
for s,f in [('TAKi','5Z'),('TAKi','PD0325901'),('MEKi','5Z'),('MEKi','PD0325901'),('PI3Ki','PI103'),('PI3Ki-TAKi','PI103+5Z'),('PI3Ki-MEKi','PI103+PD')]:
    m=pd.read_csv(f'out/tide_all_p10000_{f}.tsv',sep='\t')
    d=S[s].copy(); d['task_id']='X'+d.id.astype(str)
    j=d.merge(m,on='task_id'); print('sheet',s,'vs my',f,'n',len(j),'r=%.3f'%np.corrcoef(j.score_x,j.score_y)[0,1])
# 
d=S['PI3Ki']; print(d.sort_values('score').head(3)[['id','desc','score','p']])
for s in S:
    r=S[s][S[s].id.astype(str).isin(['138','140','141'])][['id','desc','score','p']]; print(s); print(r.to_string(index=False))
# fraction down among significant
for s in S:
    sg=S[s][S[s].p<0.025]; print(s,'sig',len(sg),'down frac %.2f'%((sg.score<0).mean()))
