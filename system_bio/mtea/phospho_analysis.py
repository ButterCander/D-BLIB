"""Phosphoproteomics (Tsirvouli repo; 0.5/2/8 h) vs the transcriptome: is enzyme regulation visible beyond mRNA?"""
import json, re
import numpy as np, pandas as pd
from scipy import stats
r = pd.read_csv("tsirvouli/All_ratios_Limma.txt"); r.columns = [c.strip() for c in r.columns]
key = r["PTM_collapse_key"]; gene = key.str.split("_").str[0]
CONDS = ["5Z", "PD", "PI", "PI5Z", "PIPD"]; TIMES = ["30m", "2h", "8h"]
LF = {(c, t): r[f"logFC_{c}_{t}"].values for c in CONDS for t in TIMES}
AP = {(c, t): r[f"adj.P.Val_{c}_{t}"].values for c in CONDS for t in TIMES}
hg = pd.read_csv("mteapy/src/mteapy/data/HumanGEM_genes.tsv", sep="\t"); met = set(hg.geneSymbol)
ismet = gene.isin(met).values
rna_name = {"5Z": "OXO", "PD": "PD", "PI": "PI", "PI5Z": "PIOXO", "PIPD": "PIPD"}
# ---------- mRNA log2FC vs DMSO at the same time (simple mean differences of log2 expression) ----------
Y = pd.read_csv("tsirvouli/counts_norm.csv", index_col=0)
gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0); sym = gi.reindex(Y.index)["gene_name"]
ok = sym.notna() & ~sym.duplicated(); Y = Y[ok.values]; Y.index = sym[ok].values
meta = pd.DataFrame([re.match(r"([A-Za-z0-9]+)_(\d+)_Rep(\d)", c).groups() for c in Y.columns], columns=["trt", "t", "rep"], index=Y.columns); meta["t"] = meta.t.astype(int)
def mrna_lfc(trt, t):
    a = Y[meta.index[(meta.trt == trt) & (meta.t == t)]].mean(axis=1); b = Y[meta.index[(meta.trt == "DMSO") & (meta.t == t)]].mean(axis=1)
    return a - b
M = {(c, t): mrna_lfc(rna_name[c], t) for c in CONDS for t in (1, 2, 8, 24)}
tmap = {"30m": 1, "2h": 2, "8h": 8}   # nearest RNA time points (1, 2, 8 h)

# ---------- (a) enzymes among significant sites ----------
print("[a] significant sites (adj.P<.05 and |logFC|>=0.5): share on metabolic enzymes vs all observed sites")
tot_share = ismet.mean()
print(f"    enzyme share among all sites: {tot_share:.3f} ({ismet.sum()} sites on {len(set(gene[ismet]))} enzymes)")
rows = []
for c in CONDS:
    row = [c]
    for t in TIMES:
        sig = (AP[(c, t)] < .05) & (np.abs(LF[(c, t)]) >= 0.5)
        obs = np.isfinite(LF[(c, t)])
        a = int((sig & ismet).sum()); b = int((sig & ~ismet).sum()); cc = int((obs & ismet & ~sig).sum()); d = int((obs & ~ismet & ~sig).sum())
        odds, p = stats.fisher_exact([[a, b], [cc, d]])
        row.append(f"{int(sig.sum()):4d} ({100*a/max(sig.sum(),1):.1f}% enz, OR {odds:.2f})")
    rows.append(row)
print("    " + f"{'cond':6}" + "".join(f"{t:>30}" for t in TIMES))
for row in rows: print("    " + f"{row[0]:6}" + "".join(f"{x:>30}" for x in row[1:]))

# ---------- (b) phospho vs mRNA on enzymes ----------
print("\n[b] enzymes with phospho data: gene-level mean site logFC vs mRNA log2FC (same condition), Spearman; and phospho-only regulation")
genes_m = sorted(set(gene[ismet]))
for c in CONDS:
    out = []
    for t in TIMES:
        s = pd.Series(LF[(c, t)], index=gene.values)[ismet]; s = s.groupby(level=0).mean().dropna()
        for rt in sorted({tmap[t], 24}):
            mm = M[(c, rt)].reindex(s.index).dropna(); ss = s.reindex(mm.index)
            out.append(f"{t}~RNA{rt}h rho={stats.spearmanr(ss, mm)[0]:+.2f}")
    print(f"    {c:5}: " + " | ".join(out))
c, t = "PIPD", "8h"
df = pd.DataFrame({"g": gene.values[ismet], "lf": LF[(c, t)][ismet], "ap": AP[(c, t)][ismet]}).dropna()
df["alf"] = df.lf.abs()
top = df.loc[df.groupby("g")["alf"].idxmax()].set_index("g")
rn = M[(c, 24)].reindex(top.index)
cand = (top.ap < .05) & (top.alf >= 0.5)
low = cand & (rn.abs() < .25)
print(f"    PIPD 8 h: enzymes with a significant phosphosite change (|logFC|>=.5): {int(cand.sum())}; of these |mRNA log2FC (24 h)|<0.25: {int(low.sum())} ({100*low.sum()/max(cand.sum(),1):.0f}%)")

# ---------- (c) task level ----------
tasks = json.load(open("out/task_genes.json"))
meta_t = pd.read_csv("mteapy/src/mteapy/data/task_metadata.tsv", sep="\t").set_index("ID")
def task_scores(vec):
    s = pd.Series(vec).groupby(level=0).mean().dropna()
    return pd.Series({t: s.reindex(g).dropna().mean() for t, g in tasks.items() if s.reindex(g).dropna().shape[0] >= 3})
print("\n[c] task-level: mean enzyme phospho logFC (8 h) vs mRNA mean log2FC (24 h), Spearman across tasks")
for c in CONDS:
    ph = task_scores(pd.Series(LF[(c, "8h")], index=gene.values)[ismet]); rn = task_scores(M[(c, 24)])
    j = pd.concat([ph, rn], axis=1, keys=["ph", "rna"]).dropna()
    print(f"    {c:5}: n tasks {len(j)}  rho = {stats.spearmanr(j.ph, j.rna)[0]:+.2f}")
ph = task_scores(pd.Series(LF[("PIPD", "8h")], index=gene.values)[ismet]); rn = task_scores(M[("PIPD", 24)])
j = pd.concat([ph, rn], axis=1, keys=["ph", "rna"]).dropna(); j["desc"] = j.index.map(meta_t.DESCRIPTION)
print("    PIPD: tasks with the largest negative phospho score (8 h) and their mRNA score (24 h):")
for i, rw in j.sort_values("ph").head(6).iterrows(): print(f"      {i:5} {rw.desc[:36]:36} phospho {rw.ph:+.2f}  mRNA {rw.rna:+.2f}")

# ---------- (d) combination vs single drugs at the phospho level ----------
print("\n[d] sites significant in the combination: share where the combination is more extreme than BOTH single drugs (same sign, margin>=0.2)")
for combo, a, b in (("PIPD", "PI", "PD"), ("PI5Z", "PI", "5Z")):
    for t in TIMES:
        lc, la, lb = LF[(combo, t)], LF[(a, t)], LF[(b, t)]
        sig = (AP[(combo, t)] < .05) & (np.abs(lc) >= 0.5) & np.isfinite(la) & np.isfinite(lb)
        beat = sig & (np.sign(lc) == np.sign(la)) & (np.sign(lc) == np.sign(lb)) & (np.abs(lc) >= np.maximum(np.abs(la), np.abs(lb)) + 0.2)
        n_m = int((sig & ismet).sum()); b_m = int((beat & ismet).sum())
        print(f"    {combo:5} {t:>4}: all sites {int(beat.sum()):4d}/{int(sig.sum()):4d} ({100*beat.sum()/max(sig.sum(),1):.0f}%) | enzymes {b_m:3d}/{n_m:3d} ({100*b_m/max(n_m,1):.0f}%)")

# ---------- (e) key enzymes ----------
print("\n[e] key enzymes: significant sites (adj.P<.05) with logFC at 8 h in PI / PD / PIPD")
for g in ["PHGDH", "PSAT1", "PSPH", "SHMT2", "ASNS", "ODC1", "AMD1", "SRM", "TXNRD1", "MTHFD1", "GLS", "ACLY"]:
    idx = np.where(gene.values == g)[0]
    if len(idx) == 0: print(f"    {g}: no phosphosite measured"); continue
    parts = []
    for i in idx:
        sig = [c for c in ("PI", "PD", "PIPD") if AP[(c, "8h")][i] < .05]
        parts.append(f"{key.iloc[i].split('_',1)[1]}[" + ", ".join(f"{c} {LF[(c,'8h')][i]:+.2f}" for c in ("PI", "PD", "PIPD")) + "]" + ("*" if sig else ""))
    print(f"    {g}: " + " ; ".join(parts[:4]))
