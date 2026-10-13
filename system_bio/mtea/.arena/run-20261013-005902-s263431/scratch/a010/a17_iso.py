from h import *
import pandas as pd
for t in ["X138","X139","X140","X141","X142"]:
    print(t, end=": ")
    for sh in auth:
        A=a_scores(sh)
        if t in A.index:
            r=A.loc[t]; print(f"{sh} {r.score:+.2f} (p {r.p:.3f}){'*' if r.p<0.025 else ''}", end=" | ")
    print(" ", a_scores("PI3Ki").loc[t,"Description"])
# HMG-CoA synthesis weak synergy common?
x = pd.ExcelFile(P("si","MOESM5.xlsx"))
for sh in ["PI3Ki-TAKi","PI3Ki-MEKi"]:
    d=x.parse(sh).iloc[:,:6]; d.columns=["task","desc","A","B","AB","dir"]
    print(sh, d[d.task.isin(["X138","X139","X140","X141","X142","X143","X144"])].round(2).values.tolist())
# MOESM7 TIDE-essential for those
e = pd.ExcelFile(P("si","MOESM7.xlsx"))
for sh in e.sheet_names[1:]:
    d=e.parse(sh); 
    s=d[d.iloc[:,0].isin(["X138","X140","X141"])]
    print(sh, s.values.tolist())
print(e.parse("Summary", header=None).to_string())
