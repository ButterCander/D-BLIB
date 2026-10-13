import pandas as pd, numpy as np
from h import *
B="D:/Cander/School/3-2/system_bio/mtea/"
x = pd.ExcelFile(B+"si/MOESM5.xlsx")
res={}
for sh,(sa,sb) in {"PI3Ki-TAKi":("PI3Ki","TAKi"),"PI3Ki-MEKi":("PI3Ki","MEKi")}.items():
    d = x.parse(sh).iloc[:,:6]; d.columns=["task","desc","A","B","AB","dir"]
    d = d.dropna(subset=["task"])
    n=len(d)
    strong = d[(d.A==0)&(d.B==0)]
    weak = d[~((d.A==0)&(d.B==0))].copy()
    # margin over best single in extremeness: for down: min(A,B)-AB ; up: AB-max(A,B)
    def margin(r):
        if r.dir.startswith("up"): return r.AB-max(r.A,r.B)
        return min(r.A,r.B)-r.AB
    weak["margin"]=weak.apply(margin,axis=1)
    print(sh,"rows",n,"strong",len(strong),"weak",len(weak))
    print("   weak margin: median %.3f  mean %.3f  frac<0.10 %.2f  frac<0.05 %.2f  n<0.1 %d"%(weak.margin.median(),weak.margin.mean(),(weak.margin<0.10).mean(),(weak.margin<0.05).mean(),(weak.margin<0.10).sum()))
    print("   strong scores:", strong[["task","desc","AB"]].round(2).values.tolist())
    # compare with authors' p of combination: join with MOESM4
    A = a_scores(sh)
    d["tid"]=d.task
    j = d.set_index("task").join(A[["p","score"]])
    print("   AB p-values (combo vs control) of calls: <.025 count", (j.p<0.025).sum(), "of", len(j), " max p", j.p.max())
    q = pd.Series(bh(np.minimum(2*A.p.values,1)), index=A.index)
    print("   calls with BH q<.05 (combo vs control):", int((q.reindex(d.task.values)<0.05).sum()), " q<.10:", int((q.reindex(d.task.values)<0.10).sum()))
    res[sh]=(d,weak,strong)
# are all synergy calls inside the significant set? "49 of 66": 
for sh in ["PI3Ki-TAKi","PI3Ki-MEKi"]:
    d=res[sh][0]; A=a_scores(sh)
    print(sh, "calls that are p<0.025 in sheet:", int((A.reindex(d.task).p<0.025).sum()), "of", len(d), "; total significant tasks", int((A.p<0.025).sum()))
# for strong: the single drug p-values
for sh,(sa,sb) in {"PI3Ki-TAKi":("PI3Ki","TAKi"),"PI3Ki-MEKi":("PI3Ki","MEKi")}.items():
    d,weak,strong=res[sh]
    for _,r in strong.iterrows():
        pa=a_scores(sa).loc[r.task,"p"]; pb=a_scores(sb).loc[r.task,"p"]; pab=a_scores(sh).loc[r.task,"p"]
        print(sh, r.task, r.desc[:40], "AB %.2f"%r.AB, "p_A %.3f p_B %.3f p_AB %.3f"%(pa,pb,pab), "score_A %.2f score_B %.2f"%(a_scores(sa).loc[r.task,"score"],a_scores(sb).loc[r.task,"score"]))
