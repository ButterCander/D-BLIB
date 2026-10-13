import openpyxl, pandas as pd, numpy as np, sys
from scipy.stats import pearsonr, spearmanr
sys.stdout.reconfigure(encoding='utf-8')
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=False)
auth={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    hdr=rows[0]
    recs=[]
    for r in rows[2:]:
        if r[0] is None: continue
        try: tid=int(r[0])
        except: continue
        recs.append(('X%d'%tid,r[1],r[2],r[3]))
    df=pd.DataFrame(recs,columns=['task','desc','score','p']).set_index('task')
    auth[ws.title]=df
    print(ws.title,'header cell:',hdr[2],'n=',len(df),'p<0.025:',(df.p<0.025).sum(),'up',((df.p<0.025)&(df.score>0)).sum(),'down',((df.p<0.025)&(df.score<0)).sum())
mine={}
for k,f in {'TAKi':'tide_5Z','MEKi':'tide_PD0325901','PI3Ki':'tide_PI103','PI3Ki-TAKi':'tide_PI103+5Z','PI3Ki-MEKi':'tide_PI103+PD'}.items():
    mine[k]=pd.read_csv(f'out/{f}.tsv',sep='\t').set_index('task_id')
# which sheet matches which of my drug?
for sh in auth:
    for k in mine:
        j=auth[sh].join(mine[k].rename(columns={'score':'score_y','pvalue':'p_y'}),how='inner').rename(columns={'score':'score_x'})
        # authors' tasks only listed when ... check n
        print(f'sheet {sh:11s} vs mine {k:11s} n={len(j)} r={pearsonr(j.score_x,j.score_y)[0]:.3f} rho={spearmanr(j.score_x,j.score_y)[0]:.3f}')
