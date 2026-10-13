import os, numpy as np, pandas as pd
from scipy.stats import false_discovery_control
OUT="out"
m=pd.read_csv("mteapy/src/mteapy/data/task_metadata.tsv",sep="\t").set_index("ID")
syn=pd.ExcelFile("si/MOESM5.xlsx"); auth={}
for combo in ("PI3Ki-TAKi","PI3Ki-MEKi"):
    d=syn.parse(combo).iloc[:,:6]; d.columns=["id","desc","A","B","AB","dir"]; d=d.dropna(subset=["id"]); d["id"]=d.id.astype(str)
    d["strong"]=(pd.to_numeric(d.A)==0)&(pd.to_numeric(d.B)==0); auth[combo]=d.set_index("id")
def load(meth,k):
    return pd.read_csv(f"{OUT}/inter_{meth}_all_p3000_{k}.tsv",sep="\t").set_index("task_id")
rows=[]
for combo,key,a,b in (("PI3Ki-MEKi","PI3Ki_MEKi","PI3Ki","MEKi"),("PI3Ki-TAKi","PI3Ki_TAKi","PI3Ki","TAKi")):
    for meth in ("tide","tidee"):
        try: d1=load(meth,f"{key}_vs{a}"); d2=load(meth,f"{key}_vs{b}")
        except FileNotFoundError: continue
        j=d1.join(d2,lsuffix="_A",rsuffix="_B",how="inner")
        j["pIU"]=np.minimum(1,2*np.maximum(j.pvalue_A,j.pvalue_B))
        j["dir"]=np.where((j.score_A<0)&(j.score_B<0),"down",np.where((j.score_A>0)&(j.score_B>0),"up","mixed"))
        j["q"]=false_discovery_control(j.pIU.values,method="bh")
        raw=j[(j.pvalue_A<.025)&(j.pvalue_B<.025)&(j.dir!="mixed")]; bh=j[(j.q<.05)&(j.dir!="mixed")]
        au=auth[combo]; ids=set(au.index); strong=set(au.index[au.strong]); weak=ids-strong
        print(f"{combo:11} {meth:5}: both contrasts p<.025 & same sign: {len(raw):3d} (down {int((raw.dir=='down').sum())}, up {int((raw.dir=='up').sum())}) | IU-BH q<.05: {len(bh):3d} | authors' calls supported raw {len(ids&set(raw.index)):2d}/{len(ids)} BH {len(ids&set(bh.index)):2d}/{len(ids)} (weak {len(weak&set(bh.index))}/{len(weak)}, strong {len(strong&set(bh.index))}/{len(strong)})")
        if meth=="tide":
            j["desc"]=j.index.map(m.DESCRIPTION)
            print("   top HSA tasks:",[(i,j.loc[i,'desc'][:26],round(j.loc[i,'score_A'],2),round(j.loc[i,'score_B'],2),round(j.loc[i,'q'],3)) for i in bh.sort_values('pIU').index[:6]])
            for k in ("X112","X17","X178","X68","X8","X76","X113","X94"):
                if k in j.index: print(f"   {k} {m.loc[k,'DESCRIPTION'][:30]:30} vs{a}: {j.loc[k,'score_A']:+.2f} (p={j.loc[k,'pvalue_A']:.3f}) vs{b}: {j.loc[k,'score_B']:+.2f} (p={j.loc[k,'pvalue_B']:.3f}) IU-q={j.loc[k,'q']:.3f}")
