"""External evidence: DepMap CRISPR dependency of AGS (24Q4) and CCLE metabolomics (Li et al. 2019)."""
import json, re
import numpy as np, pandas as pd
from scipy import stats
summ = pd.read_csv("external/crispr_gene_summary.tsv", sep="\t", index_col=0)
rows = pd.read_csv("external/crispr_gastric_rows.tsv", sep="\t", index_col=0)
model = pd.read_csv("external/Model_24Q4.csv"); agsid = "ACH-000880"
sym = lambda c: re.sub(r" \(\d+\)$", "", c)
rows.columns = [sym(c) for c in rows.columns]; summ.index = [sym(c) for c in summ.index]
ags = rows.loc[agsid]
gastric = rows.drop(index=agsid)
print(f"DepMap 24Q4: AGS {agsid}; {len(summ)} genes; {int(summ.n.max())} lines; {len(gastric)} other oesophagus/stomach lines")
z = (ags - summ["mean"]) / summ["sd"]
print(f"AGS: genes with Chronos <= -0.5: {int((ags <= -0.5).sum())}; <= -1: {int((ags <= -1).sum())}")

print("\n[1] candidate genes from the paper / my analysis: AGS Chronos (more negative = more essential), AGS z vs all lines, % of lines dependent (<-0.5), gastric-line median")
cand = ["ASNS", "ODC1", "SRM", "AMD1", "PHGDH", "PSAT1", "PSPH", "SHMT2", "MTHFD2", "TXNRD1", "COQ3", "COQ2", "COQ5", "RRM1", "RRM2", "TYMS", "GART", "PPAT", "PAICS", "ATIC", "IMPDH2", "ASS1", "GOT2", "GLS", "GLUL", "MAT2A", "DHFR", "DUT", "UMPS", "CAD", "DHODH", "ACLY", "HMGCR", "FASN", "SLC1A5", "HEXB", "GLB1", "NEU1", "GNS"]
tb = pd.DataFrame({"AGS": ags.reindex(cand), "z_vs_all": z.reindex(cand), "pct_lines_dep": 100 * summ["frac_dep"].reindex(cand), "gastric_median": gastric[cand].median()}).round(2)
print(tb.dropna().to_string())

print("\n[2] does the TIDE-essential gene concept track real essentiality? task-essential genes vs other Human-GEM genes")
hg = pd.read_csv("mteapy/src/mteapy/data/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneID").set_index("geneID")["geneSymbol"]
ess = pd.read_csv("mteapy/src/mteapy/data/HumanGEM_essential_genes_matrix.tsv", sep="\t", index_col=0).astype(bool)
ess_genes = {hg[g] for g in ess.index[ess.any(axis=1)] if g in hg.index}
allmet = {s for s in hg.values}
others = {g for g in allmet - ess_genes if g in ags.index}; ess_in = {g for g in ess_genes if g in ags.index}
for name, v in (("AGS Chronos", ags), ("pan-cancer mean Chronos", summ["mean"]), ("% lines dependent", 100 * summ["frac_dep"])):
    a = v.reindex(list(ess_in)).dropna(); b = v.reindex(list(others)).dropna()
    auc = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
    print(f"    {name:24}: task-essential genes n={len(a)} median {a.median():+.3f} | other metabolic genes n={len(b)} median {b.median():+.3f} | P(task-essential more essential)={(1-auc if 'Chronos' in name else auc):.2f}")
a_dep = (ags.reindex(list(ess_in)) <= -0.5).mean(); b_dep = (ags.reindex(list(others)) <= -0.5).mean()
print(f"    share of genes with AGS Chronos <= -0.5: task-essential {100*a_dep:.1f}% vs other metabolic {100*b_dep:.1f}%")

print("\n[3] are the authors' PI3Ki+MEKi 'synergy' tasks made of genes that matter for AGS growth? (TIDE-essential genes per task)")
syn = pd.ExcelFile("si/MOESM5.xlsx").parse("PI3Ki-MEKi").iloc[:, :2]; syn_ids = set(syn.iloc[:, 0].dropna().astype(str))
syn_ids = {("X" + i.lstrip("X")) for i in syn_ids}
def task_dep(t):
    g = [hg[x] for x in ess.index[ess[t]] if x in hg.index]; v = ags.reindex(g).dropna()
    return v.min() if len(v) else np.nan, v.mean() if len(v) else np.nan, len(v)
td = pd.DataFrame({t: task_dep(t) for t in ess.columns}, index=["min", "mean", "n"]).T.dropna()
td["syn"] = td.index.isin(syn_ids)
for col in ("min", "mean"):
    a = td.loc[td.syn, col]; b = td.loc[~td.syn, col]
    print(f"    {col:4} AGS Chronos of task-essential genes: synergy-called tasks n={len(a)} median {a.median():+.2f} | other tasks n={len(b)} median {b.median():+.2f} | MW p={stats.mannwhitneyu(a, b).pvalue:.3f}")
meta_t = pd.read_csv("mteapy/src/mteapy/data/task_metadata.tsv", sep="\t").set_index("ID")
top = td.sort_values("min").head(8)
print("    most growth-essential tasks (min Chronos):", [(i, meta_t.loc[i, "DESCRIPTION"][:26], round(r["min"], 2)) for i, r in top.iterrows()])

print("\n[4] CCLE metabolomics (Li et al. 2019): AGS vs all lines")
mb = pd.read_csv("external/CCLE_metabolomics_20190502.csv")
print(f"    {mb.shape[0]} lines x {mb.shape[1]-2} metabolites; AGS present: {'AGS_STOMACH' in set(mb.CCLE_ID)}")
a = mb[mb.CCLE_ID == "AGS_STOMACH"]
if len(a):
    a = a.iloc[0]; m = mb.drop(columns=["CCLE_ID", "DepMap_ID"])
    want = [c for c in m.columns if re.search(r"asparagine|aspartate|serine|glycine|spermidine|spermine|putrescine|ornithine|glutamine|glutamate|methionine|cysteine|proline|homoserine|S-adenosyl|NAD|dTTP|dATP|UMP|IMP|GMP|ubiquin|cholesterol|CDP", c, re.I)]
    out = []
    for c in want:
        col = m[c].astype(float); v = float(a[c])
        if np.isnan(v): continue
        pct = 100 * (col.dropna() < v).mean(); out.append((c, round(v, 2), round(pct)))
    print("    metabolite, AGS value (log2 relative), percentile among all lines:")
    for c, v, p in out: print(f"      {c[:32]:32} {v:+.2f}  {p:3d}th")
