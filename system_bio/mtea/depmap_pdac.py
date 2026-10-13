import urllib.request, time
import numpy as np, pandas as pd
t0=time.time()
model=pd.read_csv("external/Model_24Q4.csv")
pan=set(model.loc[model.OncotreeLineage.eq("Pancreas"),"ModelID"])
req=urllib.request.Request("https://ndownloader.figshare.com/files/51064667",headers={"User-Agent":"Mozilla/5.0"})
reader=pd.read_csv(urllib.request.urlopen(req,timeout=120),index_col=0,chunksize=100)
keep=[]
for ch in reader:
    keep.append(ch[ch.index.isin(pan)])
d=pd.concat(keep); d.columns=[c.split(" (")[0] for c in d.columns]
d.to_csv("external/crispr_pancreas_rows.tsv",sep="\t"); print("pancreas rows",d.shape,f"{time.time()-t0:.0f}s")
