import pandas as pd, numpy as np, warnings
from scipy.stats import false_discovery_control as bh
warnings.filterwarnings('ignore')
au={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=s,header=None).iloc[2:,[0,2,3]]; d.columns=['id','score','p']
    d['id']='X'+d.id.astype(str); d['p']=pd.to_numeric(d.p); d['score']=pd.to_numeric(d.score); au[s]=d.set_index('id')
print('AUTHORS (sheet counts p<0.025 | BH q<.05 with 2p | BH q<.10 | BH with one-tail p)')
for s,d in au.items():
    q=bh(np.clip(2*d.p.values,1e-12,1))
    q1=bh(np.clip(d.p.values,1e-12,1))
    print(s,(d.p<0.025).sum(),(q<.05).sum(),(q<.10).sum(),'| one-tail p BH<.05',(q1<.05).sum())
mine={'TAKi':'5Z','MEKi':'PD0325901','PI3Ki':'PI103','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
print('MINE (all-genes null, 10000 perm)')
for s,f in mine.items():
    m=pd.read_csv(f'out/tide_all_p10000_{f}.tsv',sep='\t')
    q=bh(np.clip(2*m.pvalue.values,1e-12,1)); print(s,(m.pvalue<0.025).sum(),(q<.05).sum(), 'n',len(m))
print('MINE (model-genes null, p10000)')
for s,f in mine.items():
    m=pd.read_csv(f'out/tide_p10000_{f}.tsv',sep='\t')
    q=bh(np.clip(2*m.pvalue.values,1e-12,1)); print(s,(m.pvalue<0.025).sum(),(q<.05).sum(),'n',len(m))
# direction share
tot=sum((d.p<0.025).sum() for d in au.values()); neg=sum(((d.p<0.025)&(d.score<0)).sum() for d in au.values())
print('neg share',neg,tot,neg/tot)
# ASNS task: Asp->Asn
for s,d in au.items():
    r=pd.read_excel('si/MOESM4.xlsx',sheet_name=s,header=None).iloc[2:,:4]
    r.columns=['id','desc','score','p']
    x=r[r.desc.astype(str).str.contains('aspartate|Asp',case=False)&r.desc.astype(str).str.contains('asparagine|Asn',case=False)]
    print(s,x[['id','desc','score','p']].values.tolist())
