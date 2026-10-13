import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
sh={}
for s in ['TAKi','MEKi']:
    d=pd.read_excel('si/MOESM4.xlsx',sheet_name=s,header=None)
    print(s,'header cell C1 =',d.iloc[0,2])
    d=d.iloc[2:,[0,2,3]]; d.columns=['id','score','p']; d['id']='X'+d.id.astype(str); d['score']=pd.to_numeric(d.score); d['p']=pd.to_numeric(d.p)
    sh[s]=d.set_index('id')
excl=['X124','X125','X126','X127','X128']
for k,f in {'5Z':'5Z','PD':'PD0325901'}.items():
    m=pd.read_csv(f'out/tide_all_p10000_{f}.tsv',sep='\t').set_index('task_id')
    for s,d in sh.items():
        j=d.join(m['score'],rsuffix='_m',how='inner'); j=j.drop(excl,errors='ignore')
        print(k,'vs sheet',s,len(j),'pearson',round(np.corrcoef(j.score,j.score_m)[0,1],3))
for t in ['X138','X140','X141']:
    print(t,{s:(round(sh[s].loc[t,'score'],3),sh[s].loc[t,'p']) for s in sh})
