import pandas as pd, numpy as np, pickle, sys
sys.stdout.reconfigure(encoding="utf-8")
from scipy.stats import false_discovery_control as bh
D = pickle.load(open(".arena/run-20261013-005902-s263431/scratch/a001/m4.pkl","rb"))
for d in D.values(): d["tid"]="X"+d.id.astype(int).astype(str)
def syn(sheet, A, B, AB):
    s = pd.read_excel("si/MOESM5.xlsx", sheet_name=sheet)
    s = s.iloc[:, :6]; s.columns=["tid","desc",A,B,AB,"dir"]
    s = s[s.tid.astype(str).str.startswith("X")]
    return s
for sheet,(A,B,AB) in {"PI3Ki-MEKi":("PI3Ki","MEKi","PI3Ki-MEKi"),"PI3Ki-TAKi":("PI3Ki","TAKi","PI3Ki-TAKi")}.items():
    s = syn(sheet,A,B,AB)
    print(sheet, "rows", len(s), "direction counts", s.dir.value_counts().to_dict())
    strong = s[(s[A]==0)&(s[B]==0)]
    print("  strong (A=B=0):", len(strong), " weak:", len(s)-len(strong))
    # combos in MOESM4
    dab = D[AB].set_index("tid"); da = D[A].set_index("tid"); db = D[B].set_index("tid")
    q = pd.Series(bh(np.clip(2*dab.p.values,1e-12,1)), index=dab.index)
    s["p_AB"]=s.tid.map(dab.p); s["q_AB"]=s.tid.map(q); s["p_A"]=s.tid.map(da.p); s["p_B"]=s.tid.map(db.p)
    s["score_AB_m4"]=s.tid.map(dab.score)
    print("  all calls: p_AB<.025:", (s.p_AB<.025).sum(), " BH q<.05:", (s.q_AB<.05).sum(), " q<.10:", (s.q_AB<.10).sum())
    st = s[(s[A]==0)&(s[B]==0)]
    print(st[["tid","desc","p_A","p_B","p_AB","q_AB"]].to_string())
    # weak margin
    w = s[~((s[A]==0)&(s[B]==0))].copy()
    sign = np.sign(w[AB])
    best = np.where(w[AB]>0, w[[A,B]].max(axis=1), w[[A,B]].min(axis=1))
    w["margin"] = (w[AB]-best).abs()
    print("  weak n", len(w), "median margin %.3f, frac<0.10: %.2f" % (w.margin.median(), (w.margin<0.10).mean()))
    # also: are single drug scores zero iff p>=.025 ?
    for nm,dd in ((A,da),(B,db)):
        z = s[s[nm]==0].tid.map(dd.p)
        print("   ", nm,"zero in synergy file -> p in MOESM4 range", z.min(), z.max(), " n=",len(z))
        nz = s[s[nm]!=0]; 
        print("   ", nm,"nonzero -> p max", nz.tid.map(dd.p).max())
