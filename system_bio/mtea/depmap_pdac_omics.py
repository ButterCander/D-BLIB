import urllib.request, time
import numpy as np, pandas as pd
t0 = time.time()
model = pd.read_csv("external/Model_24Q4.csv")
pan = set(model.loc[model.OncotreeLineage.eq("Pancreas"), "ModelID"])
def stream(fid, name):
    req = urllib.request.Request(f"https://ndownloader.figshare.com/files/{fid}", headers={"User-Agent": "Mozilla/5.0"})
    keep = []
    for ch in pd.read_csv(urllib.request.urlopen(req, timeout=180), index_col=0, chunksize=100):
        keep.append(ch[ch.index.isin(pan)])
    d = pd.concat(keep); d.columns = [c.split(" (")[0] for c in d.columns]
    d = d.loc[:, ~d.columns.duplicated()]
    d.to_csv(f"external/{name}.tsv", sep="\t"); print(name, d.shape, f"{time.time()-t0:.0f}s", flush=True)
stream(51065489, "expr_pancreas")
stream(51065747, "mut_damaging_pancreas")
stream(51065750, "mut_hotspot_pancreas")
print("done")
