import pandas as pd, numpy as np
P='si/MOESM4.xlsx'
S={}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel(P,sheet_name=s,header=None).iloc[2:,[0,2,3]]; d.columns=['id','score','p']
    d=d.dropna(subset=['p']); d['id']='X'+d.id.astype(int).astype(str); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p); S[s]=d.set_index('id')
for comb,partner in [('PI3Ki-MEKi','MEKi'),('PI3Ki-TAKi','TAKi')]:
    m=pd.read_excel('si/MOESM5.xlsx',sheet_name=comb).iloc[:,:6]; m.columns=['task_id','desc','A','B','AB','dir']
    m['AB_p']=S[comb].p.reindex(m.task_id).values; m['AB_s']=S[comb].score.reindex(m.task_id).values
    for nm,sh,col in [('PI3Ki','PI3Ki','A'),(partner,partner,'B')]:
        m[nm+'_s']=S[sh].score.reindex(m.task_id).values; m[nm+'_p']=S[sh].p.reindex(m.task_id).values
    # check masking
    for nm,col in [('PI3Ki','A'),(partner,'B')]:
        masked=np.where(m[nm+'_p']<0.025,m[nm+'_s'],0.0)
        print(comb,nm,'listed==masked score?',np.allclose(masked,m[col].astype(float),atol=1e-6),'n zero listed',(m[col].astype(float)==0).sum())
    # weak vs strong
    m['strong']=(m.A.astype(float)==0)&(m.B.astype(float)==0)
    # does the AB score exceed actual unmasked singles?
    sign=np.sign(m.AB_s)
    ext_AB=m.AB_s*sign
    ext_A=m['PI3Ki_s']*sign; ext_B=m[partner+'_s']*sign
    m['beats_both_unmasked']=(ext_AB>ext_A)&(ext_AB>ext_B)
    m['beats_PI3Ki']=ext_AB>ext_A
    m['beats_partner']=ext_AB>ext_B
    w=m[~m.strong]
    print(comb,'n calls',len(m),'strong',m.strong.sum(),'weak',len(w))
    print(' weak calls where AB beats BOTH singles using actual unmasked scores:',w.beats_both_unmasked.sum(),'of',len(w))
    print(' weak calls where one single is listed as 0 (non-sig):',((w.A.astype(float)==0)|(w.B.astype(float)==0)).sum())
    print(' strong calls beat both unmasked:',m[m.strong].beats_both_unmasked.sum(),'of',m.strong.sum())
    print(' AB significant p<0.025:',(m.AB_p<0.025).sum(),'of',len(m))
    # margin over best single, unmasked
    best=np.maximum(ext_A,ext_B); marg=(ext_AB-best)
    print(' margin over best single (unmasked, log2): median %.3f, share <0.10: %.2f, share<0: %.2f'%(marg.median(),(marg<0.10).mean(),(marg<0).mean()))
    m.to_csv(f'.arena/run-20261013-005902-s263431/scratch/a003/syn_{comb}.csv',index=False)
