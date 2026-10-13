import pandas as pd, numpy as np, pickle, sys
sys.stdout.reconfigure(encoding="utf-8")
D = pickle.load(open(".arena/run-20261013-005902-s263431/scratch/a001/m4.pkl","rb"))
for d in D.values(): d["tid"]="X"+d.id.astype(int).astype(str)
meta = pd.read_excel("si/MOESM3.xlsx", sheet_name="Task Metadata")
print(meta[meta.id.isin(["X138","X140","X141","X137","X139","X112","X17","X68","X178","X8","X13","X14","X15","X16"])][["id","SYSTEM","SUBSYSTEM","DESCRIPTION","N.GENES"]].to_string())
rows=[]
for t in ["X137","X138","X139","X140","X141","X112","X68","X178"]:
    r={"task":t, "desc":D["MEKi"].set_index("tid").loc[t,"desc"][:45]}
    for s in ["TAKi","MEKi","PI3Ki","PI3Ki-TAKi","PI3Ki-MEKi"]:
        d=D[s].set_index("tid").loc[t]; r[s]="%.2f(p%.3f)"%(d.score,d.p)
    rows.append(r)
print(pd.DataFrame(rows).to_string())
# Asp->Asn
for s in D:
    d=D[s]; x=d[d.desc.str.contains("asparagine",case=False)]
    print(s, x[["tid","desc","score","p"]].to_string(header=False))
