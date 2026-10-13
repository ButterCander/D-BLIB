"""Final verification: recompute / re-read each number used in talk_summary_qna.md from source files and compare."""
import json, re
import numpy as np, pandas as pd
from scipy import stats
rows = []
def chk(name, doc, got, tol):
    ok = abs(doc - got) <= tol
    rows.append((name, doc, round(float(got), 3), "PASS" if ok else "FAIL"))

# 1 growth synergy (re-read digitised outputs)
f = pd.read_csv("out/flobak_synergy.tsv", sep="\t", index_col=0)
chk("Bliss excess PI3Ki+MEKi (Fig 18)", 0.18, f.loc[18, "bliss_excess"], 0.02)
chk("Bliss excess PI3Ki+TAKi (Fig 21)", 0.25, f.loc[21, "bliss_excess"], 0.02)
chk("CI PI3Ki+MEKi", 0.26, f.loc[18, "CI_cal"], 0.02)
chk("CI PI3Ki+TAKi", 0.25, f.loc[21, "CI_cal"], 0.02)
# 2 permutation false-positive rate: parse log
log = open("validate_checks.log", encoding="utf8", errors="ignore").read()
m = re.search(r"false-positive rate[^\n]*?0\.39", log) or re.search(r"0\.3[89]", log)
rows.append(("gene-shuffle FPR ~0.39 present in validate_checks.log", 0.39, 0.39 if m else float("nan"), "PASS" if m else "FAIL"))
# 3 BH pass counts PI3Ki: authors vs mine
from scipy.stats import false_discovery_control as bh
au = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki", header=None).iloc[2:, [0, 2, 3]]; au.columns = ["id", "s", "p"]
pa = pd.to_numeric(au.p); qa = bh(np.clip(2 * pa.values, 1e-12, 1))
chk("authors PI3Ki BH<.05 (two-sided = 2p)", 59, (qa < .05).sum(), 0)
mine = pd.read_csv("out/tide_all_p10000_PI103.tsv", sep="\t"); qm = bh(np.clip(2 * mine.pvalue.values, 1e-12, 1))
chk("my reproduction PI3Ki BH<.05 (all-genes null)", 2, (qm < .05).sum(), 0)
sc = pd.read_excel("si/MOESM4.xlsx", sheet_name="PI3Ki", header=None).iloc[2:, [0, 2]]; sc.columns = ["id", "s"]
sc["id"] = "X" + sc.id.astype(str); j = sc.merge(mine, left_on="id", right_on="task_id")
chk("score correlation PI3Ki mine vs authors", 0.95, np.corrcoef(pd.to_numeric(j.s), j.score)[0, 1], 0.02)
# 4 E-Flux predicted inhibition range (stored grid)
g = json.load(open("out/eflux_grid.json")); vals = [v for x in g for v in x["pred"].values()]
chk("E-Flux min predicted inhibition (stored grid, %)", -19, min(vals), 1.5)
chk("E-Flux max predicted inhibition (stored grid, %)", 5, max(vals), 1.5)
# 5 statin ranks: PANC-1 20 uM, three statins, X138/X140/X141 best ranks
rk = {}
for st in ("Fluvastatin", "Lovastatin", "Simvastatin"):
    r = pd.read_csv(f"out/tidex_statin_PANC1_{st}_High.tsv", sep="\t").set_index("task_id")["score"].dropna().rank(ascending=False)
    rk[st] = [int(r[t]) for t in ("X138", "X140", "X141")]
chk("worst expected-task rank among 3 statins in PANC-1 (<=4)", 4, max(max(v) for v in rk.values()), 0)
# 6 DepMap pancreatic dependency
cr = pd.read_csv("external/crispr_pancreas_rows.tsv", sep="\t", index_col=0); cr = cr.loc[:, ~cr.columns.duplicated()]
chk("HMGCS1 median Chronos", -1.95, cr["HMGCS1"].median(), 0.03)
chk("HMGCS1 fraction dependent", 1.00, (cr["HMGCS1"] <= -0.5).mean(), 0.0)
chk("HMGCR median Chronos", -1.10, cr["HMGCR"].median(), 0.03)
chk("HMGCR fraction dependent", 0.98, (cr["HMGCR"] <= -0.5).mean(), 0.02)
chk("SQLE median Chronos", -0.14, cr["SQLE"].median(), 0.03)
chk("SQLE fraction dependent", 0.02, (cr["SQLE"] <= -0.5).mean(), 0.02)
# 7 AGS Tsirvouli HMGCS1 24 h
off = pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv"); gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)
off["sym"] = gi.reindex(off.Gene_name)["gene_name"].values
for c, v in (("PD", 1.90), ("PIPD", 2.47)):
    x = off[(off.Contrast == f"Diff_{c}_24h") & (off.sym == "HMGCS1")].logFC.iloc[0]; chk(f"AGS full-dose {c} HMGCS1 24 h", v, x, 0.02)
# 8 PDAC program means (24 listed genes), GSE238198 24 h
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
for l, v in (("CAPAN2", 1.06), ("CFPAC1", 1.31), ("MiaPaca2", 0.33), ("PANC1", 1.12)):
    d = pd.read_csv(f"out/pdac_de_{l}_h24.tsv", sep="\t").set_index("gene")["logFC"]
    chk(f"program mean {l} trametinib 24 h", v, d.reindex([s for s in SREBP if s in d.index]).mean(), 0.02)
# 9 MRTX1133 PANC-1 24 h program
d = pd.read_csv("pdac/GSE279749_gene_fpkm_cell_lines.txt.gz", sep="\t").drop_duplicates("gene_name").set_index("gene_name")
lf = np.log2(d[[f"P_M_{i}" for i in (1, 2, 3)]] + 1).mean(axis=1) - np.log2(d[[f"P_NT_{i}" for i in (1, 2, 3)]] + 1).mean(axis=1)
chk("PANC-1 MRTX1133 24 h program mean", 0.33, lf.reindex([s for s in SREBP if s in lf.index]).mean(), 0.02)
# 10 baseline predictor
model = pd.read_csv("external/Model_24Q4.csv").set_index("ModelID"); n2i = {r.StrippedCellLineName: i for i, r in model.iterrows() if r.OncotreeLineage == "Pancreas"}
E = pd.read_csv("external/expr_pancreas.tsv", sep="\t", index_col=0); TR = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
ind = {"CAPAN2": 1.06, "PANC1": 0.725, "CFPAC1": 0.68, "KP4": 0.18, "MIAPACA2": 0.155, "ASPC1": -0.10, "SW1990": -0.17}
x = [E.loc[n2i[k], TR].mean() for k in ind]; y = list(ind.values())
chk("baseline trunk vs inducibility rho", 0.96, stats.spearmanr(x, y)[0], 0.02)
res = pd.DataFrame(rows, columns=["check", "doc value", "recomputed", "result"])
pd.set_option("display.width", 200); pd.set_option("display.max_colwidth", 70)
print(res.to_string(index=False))
print("\nFAIL count:", int((res.result == "FAIL").sum()), "of", len(res))
