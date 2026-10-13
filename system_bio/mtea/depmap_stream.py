"""Stream DepMap 24Q4 CRISPRGeneEffect.csv (figshare mirror of the public release): keep AGS + gastric rows and per-gene summaries."""
import io, sys, time, urllib.request
import numpy as np, pandas as pd
t0 = time.time()
model = pd.read_csv("https://ndownloader.figshare.com/files/51065297")
model.to_csv("external/Model_24Q4.csv", index=False)
print("Model rows", len(model), list(model.columns)[:12], flush=True)
name_col = "StrippedCellLineName" if "StrippedCellLineName" in model else "CellLineName"
ags = model[model[name_col].astype(str).str.upper() == "AGS"]
print("AGS rows:", ags[["ModelID", "CellLineName", "OncotreeLineage", "OncotreePrimaryDisease"]].to_dict("records"), flush=True)
agsid = ags.ModelID.iloc[0]
gastric = set(model.loc[model.OncotreeLineage.astype(str).str.contains("Esophagus/Stomach"), "ModelID"])
n = 0; s = None; ss = None; low = None; keep = {}
req = urllib.request.Request("https://ndownloader.figshare.com/files/51064667", headers={"User-Agent": "Mozilla/5.0"})
resp = urllib.request.urlopen(req, timeout=120)
reader = pd.read_csv(resp, index_col=0, chunksize=100)
cols = None
for ch in reader:
    cols = ch.columns
    v = ch.values.astype(float)
    if s is None: s = np.zeros(v.shape[1]); ss = np.zeros(v.shape[1]); low = np.zeros(v.shape[1]); cnt = np.zeros(v.shape[1])
    ok = ~np.isnan(v)
    s += np.where(ok, v, 0).sum(0); ss += np.where(ok, v, 0).__pow__(2).sum(0); low += ((v < -0.5) & ok).sum(0); cnt += ok.sum(0)
    for mid in ch.index:
        if mid == agsid or mid in gastric: keep[mid] = ch.loc[mid]
    n += len(ch)
    if n % 400 == 0: print(f"  rows {n} ({time.time()-t0:.0f}s)", flush=True)
mean = s / cnt; sd = np.sqrt(np.maximum(ss / cnt - mean ** 2, 0))
summ = pd.DataFrame({"mean": mean, "sd": sd, "frac_dep": low / cnt, "n": cnt}, index=cols)
summ.to_csv("external/crispr_gene_summary.tsv", sep="\t")
pd.DataFrame(keep).T.to_csv("external/crispr_gastric_rows.tsv", sep="\t")
print("done rows", n, "gastric kept", len(keep), "AGS in kept:", agsid in keep, f"{time.time()-t0:.0f}s", flush=True)
