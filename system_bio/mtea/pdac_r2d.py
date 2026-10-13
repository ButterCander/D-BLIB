"""Response-to-dependency (R2D) map for pancreatic cancer.
Idea: a targeted drug (here the MEK inhibitor trametinib) changes the expression of metabolic genes; genes pushed UP that the cancer also needs
(DepMap CRISPR dependency in pancreatic lines) are candidate 'escape routes' for a partner drug; genes pushed DOWN that are essential are
candidate 'amplified vulnerabilities'.  Data: GSE238198 (4 PDAC lines, trametinib 24 h), DepMap 24Q4 pancreas (47 lines), AGS MEKi (GSE285616,
Tsirvouli 24 h), PDAC statin data (GSE149566; PANC-1, MiaPaCa-2).  Hypothesis generation only."""
import json
import numpy as np, pandas as pd
from scipy import stats
D = "mteapy/src/mteapy/data"
hg = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol")
met = set(hg.geneSymbol); id2sym = hg.set_index("geneID")["geneSymbol"]
LINES = ["CAPAN2", "CFPAC1", "MiaPaca2", "PANC1"]
# ---- response in PDAC (trametinib 24 h)
lf = {}; ap = {}
for l in LINES:
    d = pd.read_csv(f"out/pdac_de_{l}_h24.tsv", sep="\t").set_index("gene")
    lf[l] = d["logFC"]; ap[l] = d["adj.P.Val"]
LF = pd.DataFrame(lf); AP = pd.DataFrame(ap)
R = pd.DataFrame({"mean_lfc": LF.mean(axis=1), "n_lines": LF.notna().sum(axis=1)})
R["n_up"] = ((LF >= 0.5) & (AP < 0.05)).sum(axis=1); R["n_down"] = ((LF <= -0.5) & (AP < 0.05)).sum(axis=1)
R = R[R.index.isin(met) & (R.n_lines >= 3)]
# ---- dependency in pancreatic lines (Chronos), vs pan-cancer
cr = pd.read_csv("external/crispr_pancreas_rows.tsv", sep="\t", index_col=0); cr = cr.loc[:, ~cr.columns.duplicated()]
summ = pd.read_csv("external/crispr_gene_summary.tsv", sep="\t", index_col=0); summ.index = [c.split(" (")[0] for c in summ.index]
dep = pd.DataFrame({"pdac_median": cr.median(), "pdac_frac_dep": (cr <= -0.5).mean(), "pan_mean": summ["mean"].reindex(cr.columns)})
dep["selectivity"] = dep.pdac_median - dep.pan_mean       # negative = more essential in pancreas than average
T = R.join(dep, how="inner")
# ---- cross-cancer and mechanism checks
ags_g = pd.read_csv("out/de_PD0325901.tsv", sep="\t", index_col=0)["log2FoldChange"]
off = pd.read_csv("tsirvouli_official/DEGs_perContrast_all.csv"); gi = pd.read_csv("tsirvouli/gene_info.tsv", sep="\t", index_col=0)
off = off[off.Contrast == "Diff_PD_24h"].copy(); off["sym"] = gi.reindex(off.Gene_name)["gene_name"].values
ags_t = off.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")["logFC"]
T["ags_gse"] = ags_g.reindex(T.index); T["ags_tsir"] = ags_t.reindex(T.index)
st = {}
for l, nm in (("PANC1", "PANC1_Simvastatin_High"), ("MiaPaca2", "MiaPaca2_Simvastatin_High")):
    d = pd.read_csv(f"out/statin_de_{nm}.tsv", sep="\t"); s = pd.read_csv("pdac/statin_symbols.tsv", sep="\t").set_index("gene_id")["gene_short_name"]
    d["sym"] = d.gene.map(s); st[l] = d.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym")["logFC"]
T["statin"] = pd.DataFrame(st).mean(axis=1).reindex(T.index)
# ---- classes
up = (T.mean_lfc >= 0.5) & (T.n_up >= 3)
down = (T.mean_lfc <= -0.5) & (T.n_down >= 3)
dep_ok = (T.pdac_median <= -0.5) & (T.pdac_frac_dep >= 0.5)
T["class"] = np.where(up & dep_ok, "A_escape (up & essential)", np.where(down & dep_ok, "B_amplified (down & essential)", np.where(up, "up_nonessential", np.where(down, "down_nonessential", ""))))
print(f"metabolic genes with response and dependency data: {len(T)}")
print(T["class"].value_counts().to_string())
DRUG = {"HMGCR": "statins", "FDPS": "N-bisphosphonates (zoledronate)", "FNTA": "farnesyltransferase inhibitors", "FNTB": "farnesyltransferase inhibitors", "SQLE": "SQLE inhibitors", "HMGCS1": "-",
        "DHODH": "brequinar / teriflunomide", "RRM1": "gemcitabine, hydroxyurea", "RRM2": "gemcitabine, hydroxyurea", "TYMS": "5-FU, pemetrexed", "DHFR": "methotrexate, pemetrexed",
        "IMPDH2": "mycophenolate", "IMPDH1": "mycophenolate", "CAD": "PALA", "ODC1": "DFMO", "PHGDH": "NCT-503 (tool)", "GLS": "telaglenastat (CB-839)", "NAMPT": "FK866",
        "ACLY": "bempedoic acid", "FASN": "TVB-2640", "SCD": "SCD inhibitors", "TXNRD1": "auranofin", "SLC7A11": "erastin / sulfasalazine", "GPX4": "RSL3 (tool)", "LDHA": "LDHA inhibitors",
        "HK2": "2-DG / lonidamine", "MTHFD2": "MTHFD2 inhibitors", "CTPS1": "CTPS1 inhibitors", "PRPS1": "-", "GART": "lometrexol", "ATIC": "-", "PPAT": "-", "UMPS": "-", "TK1": "-"}
T["drug_class"] = [DRUG.get(g, "") for g in T.index]
cols = ["mean_lfc", "n_up", "n_down", "pdac_median", "pdac_frac_dep", "selectivity", "ags_gse", "ags_tsir", "statin", "drug_class"]
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 34)
A = T[T["class"].str.startswith("A")].sort_values("mean_lfc", ascending=False)
B = T[T["class"].str.startswith("B")].sort_values("mean_lfc")
print("\nClass A (up, essential): n =", len(A)); print(A[cols].head(25).round(2).to_string())
print("\nClass B (down, essential): n =", len(B)); print(B[cols].head(25).round(2).to_string())
# concordance
for nm, X in (("A", A), ("B", B)):
    sg = np.sign(X.mean_lfc)
    print(f"class {nm}: same direction in AGS MEKi (GSE285616) {(np.sign(X.ags_gse)==sg).mean():.2f}, in Tsirvouli {(np.sign(X.ags_tsir)==sg).mean():.2f}; "
          f"statin PDAC logFC>0.5 in {(X.statin>=0.5).mean():.2f} of genes (all metabolic genes: {(T.statin>=0.5).mean():.2f})")
iso = ["HMGCR", "HMGCS1", "MVK", "PMVK", "MVD", "IDI1", "FDPS", "FDFT1", "SQLE", "LSS", "CYP51A1", "DHCR7", "DHCR24", "SC5D", "LDLR", "INSIG1", "ACAT2", "NSDHL", "MSMO1", "HSD17B7", "EBP", "TM7SF2", "SREBF2"]
print("\nIsoprenoid / sterol genes (trametinib 24 h mean logFC; PDAC median Chronos; simvastatin logFC):")
print(T.reindex([g for g in iso if g in T.index])[["mean_lfc", "n_up", "pdac_median", "pdac_frac_dep", "statin", "ags_gse"]].round(2).to_string())
T.to_csv("out/pdac_r2d.tsv", sep="\t")
