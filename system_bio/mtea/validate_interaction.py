"""Re-test the interaction / HSA contrasts with the calibrated CAMERA-like set test (estimated inter-gene correlation)."""
import json
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control
OUT = "out"; D = "mteapy/src/mteapy/data"
cnt = pd.read_csv("GSE285616_counts.csv.gz", index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]; cnt = cnt.loc[cnt.sum(axis=1) >= 10]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
Y = np.log2(x / sf + 1)
batch = np.array([i // 6 for i in range(24)]); cond = np.array([i % 6 for i in range(24)])
X = np.column_stack([np.ones(24)] + [(batch == b).astype(float) for b in (1, 2, 3)] + [(cond == c).astype(float) for c in range(1, 6)])
beta = np.linalg.lstsq(X, Y.T, rcond=None)[0]; R = (Y.T - X @ beta)
Rz = (R - R.mean(0)).T; Rz /= (np.linalg.norm(Rz, axis=1, keepdims=True) + 1e-12)
genes = list(cnt.index); gi = {g: i for i, g in enumerate(genes)}
tasks = json.load(open(f"{OUT}/task_genes.json"))
meta = pd.read_csv(f"{D}/task_metadata.tsv", sep="\t").set_index("ID")
sets = {t: np.array(sorted({gi[g] for g in gl if g in gi})) for t, gl in tasks.items()}
sets = {t: v for t, v in sets.items() if 5 <= len(v) <= 200}
rho = {t: max(0.0, (lambda C, m: (C.sum() - m) / (m * (m - 1)))(Rz[v] @ Rz[v].T, len(v))) for t, v in sets.items()}
syn = pd.ExcelFile("si/MOESM5.xlsx"); auth = {}
for combo in ("PI3Ki-TAKi", "PI3Ki-MEKi"):
    d = syn.parse(combo).iloc[:, :6]; d.columns = ["id", "desc", "A", "B", "AB", "dir"]; auth[combo] = set(d.dropna(subset=["id"]).id.astype(str))


def test(file):
    r = pd.read_csv(f"{OUT}/{file}.tsv", sep="\t", index_col=0).reindex(genes)
    z = r["stat"].values; ok = np.isfinite(z); z = np.where(ok, z, 0.0)
    rows = []
    for t, idx in sets.items():
        inset = np.zeros(len(z), bool); inset[idx] = True; m = len(idx)
        zin = z[idx][ok[idx]].mean(); zout = z[(~inset) & ok]
        vif = 1 + (m - 1) * rho[t]
        tt = (zin - zout.mean()) / np.sqrt(zout.var(ddof=1) * (vif / m + 1 / len(zout)))
        rows.append((t, zin, tt, 2 * stats.t.sf(abs(tt), len(zout) - 2)))
    d = pd.DataFrame(rows, columns=["task", "mean_z", "t", "p"]).set_index("task"); d["q"] = false_discovery_control(d.p.values, method="bh")
    return d


print("CAMERA-like re-test of the contrasts (120 task gene sets with 5-200 genes, Wald z of DESeq2, estimated rho)")
for combo, key, a in (("PI3Ki-MEKi", "PI3Ki_MEKi", "MEKi"), ("PI3Ki-TAKi", "PI3Ki_TAKi", "TAKi")):
    for f, lab in ((f"inter_{key}_dose", "dose-scaled interaction AB-.5A-.5B"), (f"inter_{key}_add", "additive interaction AB-A-B"),
                   (f"inter_{key}_vsPI3Ki", "AB vs PI3Ki alone"), (f"inter_{key}_vs{a}", f"AB vs {a} alone")):
        d = test(f)
        neg = int(((d.p < .05) & (d.mean_z < 0)).sum()); pos_ = int(((d.p < .05) & (d.mean_z > 0)).sum())
        ids = {"X" + str(i) for i in range(0, 300)}
        au = {("X" + i.lstrip("X")) for i in auth[combo]}
        hit = set(d.index[d.p < .05]) & au; hitq = set(d.index[d.q < .05]) & au
        print(f"  {combo:11} {lab:36} p<.05: {int((d.p < .05).sum()):3d} (neg {neg}, pos {pos_}) | BH q<.05: {int((d.q < .05).sum()):3d} | authors' synergy ids among p<.05: {len(hit)}, among BH: {len(hitq)} (of {len(au & set(d.index))} testable)")
        if f.endswith("_dose") and combo == "PI3Ki-MEKi":
            top = d.sort_values("p").head(8)
            print("     top:", [(i, meta.loc[int(str(i).lstrip("X")) if False else i, "DESCRIPTION"][:24], round(r.mean_z, 2), round(r.q, 3)) for i, r in top.iterrows()])
