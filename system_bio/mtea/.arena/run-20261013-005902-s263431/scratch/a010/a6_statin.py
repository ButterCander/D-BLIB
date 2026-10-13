import pandas as pd, numpy as np, glob, os
B="D:/Cander/School/3-2/system_bio/mtea/"
md = pd.read_excel(B+"si/MOESM3.xlsx", sheet_name="Task Metadata")
md = md.set_index('id')
for t in ["X138","X139","X140","X141","X76","X8","X1"]:
    print(t, md.loc[t,'SYSTEM'],'|',md.loc[t,'SUBSYSTEM'],'|',md.loc[t,'DESCRIPTION'],'| genes',md.loc[t,'N.GENES'])
print()
for f in sorted(glob.glob(B+"out/tidex_statin_*.tsv")):
    d = pd.read_csv(f, sep="\t").set_index("task_id")
    print(os.path.basename(f), d.columns.tolist()) if "PANC1_Fluvastatin_High" in f else None
    sc = d["score"].dropna()
    r = sc.rank(ascending=False, method="min")
    n = len(sc)
    out = []
    for t in ["X138","X140","X141"]:
        p = d.loc[t,"pvalue"]
        out.append(f"{t}: rank {int(r[t]):3d} score {sc[t]:+.2f} p {p:.3f}")
    print(os.path.basename(f).replace("tidex_statin_","").replace(".tsv",""), "n=",n, "|", " | ".join(out))
