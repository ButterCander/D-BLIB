import glob,os,numpy as np,pandas as pd
EXP=["X138","X140","X141"]
for f in sorted(glob.glob("out/tidex_statin_*.tsv")):
    n=os.path.basename(f)[len("tidex_statin_"):-4]
    r=pd.read_csv(f,sep="\t").set_index("task_id"); s=r.score.dropna(); rk=s.rank(ascending=False,method='min')
    print(f"{n:30} n={len(s)} ", " ".join(f"{t}: score {r.score.get(t):+.2f} rank {int(rk[t])} p {r.pvalue.get(t):.3f}" for t in EXP))
