import pandas as pd, numpy as np
from scipy.stats import false_discovery_control as bh, spearmanr, pearsonr
m={'PI3Ki':'PI103','MEKi':'PD0325901','TAKi':'5Z','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
def auth(sh):
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=sh,header=None).iloc[2:,:4]; d.columns=['id','desc','score','p']
    d=d.dropna(subset=['id']); d['id']='X'+d.id.astype(int).astype(str); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p); return d.set_index('id')
print('cond | all-genes null p<.025 -> BH | model-genes null | r(scores) vs authors | TIDE-e all-genes')
for sh,c in m.items():
    a=auth(sh)
    out=[]
    for f in ['tide_all_p10000','tide_p10000']:
        d=pd.read_csv(f'out/{f}_{c}.tsv',sep='\t').set_index('task_id')
        q=bh(np.clip(2*d.pvalue.values,1e-12,1))
        out.append('%d -> %d'%((d.pvalue<.025).sum(),(q<.05).sum()))
    d=pd.read_csv(f'out/tide_all_p10000_{c}.tsv',sep='\t').set_index('task_id')
    j=a.join(d[['score']],rsuffix='_mine',how='inner')
    r=pearsonr(j.score,j.score_mine)[0]; sp=spearmanr(j.score,j.score_mine)[0]
    e=pd.read_csv(f'out/tidee_all_p10000_{c}.tsv',sep='\t').set_index('task_id')
    qe=bh(np.clip(2*e.pvalue.values,1e-12,1))
    print('%-11s'%sh,out,'r=%.2f rho=%.2f n=%d'%(r,sp,len(j)),'TIDE-e %d -> %d (n=%d)'%((e.pvalue<.025).sum(),(qe<.05).sum(),len(e)))
