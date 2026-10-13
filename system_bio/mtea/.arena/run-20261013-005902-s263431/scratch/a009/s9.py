import openpyxl, pandas as pd, numpy as np, sys, itertools
from scipy.stats import pearsonr, spearmanr
sys.stdout.reconfigure(encoding='utf-8')
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True)
auth={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True)); recs=[]
    for r in rows[2:]:
        try: recs.append(('X%d'%int(r[0]),r[2],r[3]))
        except: pass
    auth[ws.title]=pd.DataFrame(recs,columns=['task','score','p']).set_index('task')
fm={'TAKi':'5Z','MEKi':'PD0325901','PI3Ki':'PI103','PI3Ki-TAKi':'PI103+5Z','PI3Ki-MEKi':'PI103+PD'}
mine={k:pd.read_csv(f'out/tide_{v}.tsv',sep='\t').set_index('task_id') for k,v in fm.items()}
A=pd.DataFrame({k:auth[k].score for k in fm}); M=pd.DataFrame({k:mine[k].score for k in fm}).reindex(A.index)
print('all-task r per condition:',{k:round(pearsonr(A[k],M[k])[0],3) for k in fm})
print('all 187x5 pooled r',round(pearsonr(A.values.ravel(),M.values.ravel())[0],3),'rho',round(spearmanr(A.values.ravel(),M.values.ravel())[0],3))
# candidate: six tasks
for tasks in [['X76','X113','X112','X178','X17','X8'],['X76','X113','X112','X178','X17','X68'],['X76','X113','X17','X112','X178','X68','X8']]:
    a=A.loc[tasks].values.ravel(); m=M.loc[tasks].values.ravel()
    print(tasks,'n',len(a),'r',round(pearsonr(a,m)[0],4))
# top-N tasks by |authors' PI3Ki score|
for n in [6,10,20]:
    t=A['PI3Ki'].abs().sort_values(ascending=False).index[:n]
    a=A.loc[t].values.ravel(); m=M.loc[t].values.ravel(); print('top',n,'by |PI3Ki|',list(t),'r',round(pearsonr(a,m)[0],4))
# pooled r after removing the ASNS row
t=[x for x in A.index if x!='X76']
a=A.loc[t].values.ravel(); m=M.loc[t].values.ravel(); print('pooled w/o X76 r',round(pearsonr(a,m)[0],3))
print(A.loc['X76'].round(2).to_dict(), M.loc['X76'].round(2).to_dict())
