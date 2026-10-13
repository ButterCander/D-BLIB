"""What distinguishes PDAC lines that induce the mevalonate trunk after MAPK-pathway inhibition? (n = 7 lines; exploratory)"""
import numpy as np, pandas as pd
from scipy import stats
model = pd.read_csv("external/Model_24Q4.csv").set_index("ModelID")
name2id = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
LINES = {"CAPAN2": "CAPAN2", "CFPAC1": "CFPAC1", "PANC1": "PANC1", "MIAPACA2": "MIAPACA2", "KP4": "KP4", "ASPC1": "ASPC1", "SW1990": "SW1990"}
# program mean log2FC (24 sterol/SREBP genes) after MAPK-pathway inhibition, from sections 18-19
IND = {"CAPAN2": [1.06], "CFPAC1": [1.31, 0.05], "PANC1": [1.12, 0.33], "MIAPACA2": [0.33, -0.02], "KP4": [0.18], "ASPC1": [-0.10], "SW1990": [-0.17]}
ind = pd.Series({k: np.mean(v) for k, v in IND.items()}).sort_values(ascending=False)
ids = {k: name2id[k] for k in LINES}
print("inducibility (mean program log2FC over MAPK-pathway inhibitor datasets):"); print(ind.round(2).to_string())
E = pd.read_csv("external/expr_pancreas.tsv", sep="\t", index_col=0)
M = pd.read_csv("external/mut_damaging_pancreas.tsv", sep="\t", index_col=0)
H = pd.read_csv("external/mut_hotspot_pancreas.tsv", sep="\t", index_col=0)
C = pd.read_csv("external/crispr_pancreas_rows.tsv", sep="\t", index_col=0); C = C.loc[:, ~C.columns.duplicated()]
sel = [ids[k] for k in ind.index]
print("\n[1] genotype of the 7 lines (damaging mutation = 1; KRAS hotspot allele)")
rows = []
for k in ind.index:
    i = ids[k]
    kras = [c for c in H.columns if c.startswith("KRAS") and H.loc[i, c] > 0] if i in H.index else []
    rows.append((k, round(ind[k], 2), ",".join(kras) or "-", *[int(M.loc[i, g]) if (i in M.index and g in M.columns) else -1 for g in ("TP53", "SMAD4", "CDKN2A", "ARID1A", "KRAS")]))
print(pd.DataFrame(rows, columns=["line", "induction", "KRAS hotspot", "TP53", "SMAD4", "CDKN2A", "ARID1A", "KRAS dam"]).to_string(index=False))
# a priori panels
panels = {
 "SREBP machinery": ["SREBF1", "SREBF2", "SCAP", "INSIG1", "INSIG2", "MBTPS1", "MBTPS2", "LDLR", "PCSK9", "HMGCR", "HMGCS1", "SQLE", "MVK", "MVD", "FDPS", "SCD"],
 "MAPK output / pathway": ["DUSP4", "DUSP6", "ETV4", "ETV5", "SPRY2", "FOSL1", "MYC", "KRAS", "BRAF", "NF1", "EGFR", "ERBB2", "ERBB3"],
 "subtype / EMT": ["GATA6", "HNF1A", "HNF4A", "TFF1", "AGR2", "LGALS4", "CLDN18", "KRT5", "KRT6A", "S100A2", "TP63", "LAMC2", "VIM", "ZEB1", "CDH1"],
}
print("\n[2] baseline expression (log2 TPM+1, DepMap) vs inducibility, Spearman across 7 lines (|rho| >= 0.71 flagged; exact p for n=7: 0.76 -> .05)")
out = []
for pn, gl in panels.items():
    for g in gl:
        if g in E.columns:
            v = E.loc[sel, g].values; r, p = stats.spearmanr(ind.values, v); out.append((pn, g, round(r, 2), round(p, 3), round(float(np.nanmedian(E[g])), 2)))
df = pd.DataFrame(out, columns=["panel", "gene", "rho", "p", "median in 55 pancreatic models"])
pd.set_option("display.width", 200)
print(df.reindex(df.rho.abs().sort_values(ascending=False).index).head(14).to_string(index=False))
# classical vs basal-like score
cl = ["GATA6", "HNF1A", "HNF4A", "TFF1", "AGR2", "LGALS4", "CLDN18"]; ba = ["KRT5", "KRT6A", "S100A2", "TP63", "LAMC2"]
z = (E - E.mean()) / E.std()
sc = z[[g for g in cl if g in z.columns]].mean(axis=1) - z[[g for g in ba if g in z.columns]].mean(axis=1)
print("\n[3] classical-minus-basal score (z over pancreatic models):", {k: round(float(sc[ids[k]]), 2) for k in ind.index}, "| rho =", round(stats.spearmanr(ind.values, sc[sel].values)[0], 2))
trunk = ["HMGCS1", "HMGCR", "MVK", "MVD", "IDI1", "FDPS"]
bt = E.loc[sel, trunk].mean(axis=1)
print("[3b] baseline trunk expression:", {k: round(float(v), 2) for k, v in zip(ind.index, bt.values)}, "| rho =", round(stats.spearmanr(ind.values, bt.values)[0], 2))
print("\n[4] CRISPR dependency (Chronos) vs inducibility, a priori genes")
dep = []
for g in ["KRAS", "RAF1", "BRAF", "MAP2K1", "MAPK1", "MAPK3", "SREBF2", "SCAP", "MBTPS1", "MBTPS2", "HMGCR", "HMGCS1", "MVK", "FNTA", "GGPS1", "PIK3CA", "AKT1", "MTOR", "MYC", "SHOC2", "SOS1"]:
    if g in C.columns:
        ok_ = [i for i in sel if i in C.index]; v = C.loc[ok_, g].values; r, p = stats.spearmanr(ind.reindex([k for k in ind.index if ids[k] in C.index]).values, v); dep.append((g, round(r, 2), round(p, 3), round(float(np.nanmedian(C[g])), 2)))
print(pd.DataFrame(dep, columns=["gene", "rho", "p", "median Chronos"]).to_string(index=False))

# ---------- robustness of the baseline-trunk predictor ----------
print("\n[5] empirical check of baseline trunk mean vs inducibility")
rng = np.random.default_rng(3)
pool = [g for g in E.columns if E[g].median() > 2]
obs = stats.spearmanr(ind.values, bt.values)[0]
rs = np.array([stats.spearmanr(ind.values, E.loc[sel, list(rng.choice(pool, 6, replace=False))].mean(axis=1).values)[0] for _ in range(3000)])
print(f"  observed rho {obs:+.2f}; random 6-gene sets: fraction with rho >= observed {(rs >= obs).mean():.4f}; fraction with |rho| >= |observed| {(np.abs(rs) >= abs(obs)).mean():.4f}; 95th pct of |rho| {np.percentile(np.abs(rs), 95):.2f}")
loo = []
for k in ind.index:
    keep = [x for x in ind.index if x != k]
    loo.append(round(stats.spearmanr(ind[keep].values, bt.reindex([ids[x] for x in keep]).values)[0], 2))
print("  leave-one-line-out rho:", dict(zip(ind.index, loo)))
# same relation inside the original datasets (baseline = DMSO/NT of that dataset)
print("\n[6] baseline trunk (mean log2 expression of the 6 trunk genes in the control samples of the same study) vs induction inside each study")
TR = ["HMGCS1", "HMGCR", "MVK", "MVD", "IDI1", "FDPS"]
cnt = {}
for l in ("CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"):
    d = pd.read_csv(f"pdac/GSE238198_Pdac.{l}.genes.normalizedCounts.txt.gz", sep="	", index_col=0)
    cnt[l] = np.log2(d[[c for c in d.columns if "_NT_" in c]] + 1).mean(axis=1).reindex(TR).mean()
ind238 = {"CAPAN2": 1.06, "CFPAC1": 1.31, "MiaPaca2": 0.33, "PANC1": 1.12}
print("  GSE238198 baseline:", {k: round(v, 2) for k, v in cnt.items()}, "| rho vs 24 h induction:", round(stats.spearmanr([cnt[k] for k in ind238], [ind238[k] for k in ind238])[0], 2))
d = pd.read_csv("pdac/GSE279749_gene_fpkm_cell_lines.txt.gz", sep="	").drop_duplicates("gene_name").set_index("gene_name")
b279 = {n: np.log2(d[[f"{k}_NT_{i}" for i in (1, 2, 3)]] + 1).mean(axis=1).reindex(TR).mean() for n, k in (("AsPC-1", "A"), ("PANC-1", "P"), ("SW1990", "S9"), ("CFPAC-1", "C"))}
i279 = {"AsPC-1": -0.10, "PANC-1": 0.33, "SW1990": -0.17, "CFPAC-1": 0.05}
print("  GSE279749 baseline:", {k: round(v, 2) for k, v in b279.items()}, "| rho vs 24 h induction:", round(stats.spearmanr([b279[k] for k in i279], [i279[k] for k in i279])[0], 2))
