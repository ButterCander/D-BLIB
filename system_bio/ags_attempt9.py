"""E19 screening efficiency (synergies found in the top k); E20 does the consensus allow fewer evolution runs (50 instead of 150)?"""
import warnings
from scipy.stats import rankdata
from ags_eval import *

warnings.filterwarnings("ignore")
FOLD = ["link-only", "topology-only", "topo-and-link"]
rk = rankdata
exn = lambda f, k, n: load(f, k, n, "bliss")["Response excess over subset"].reindex(pairs).fillna(0).values
vfn = lambda f, k, n: frac(load(f, k, n, "bliss", "modelwise")).fillna(0).values
auc = lambda s: roc_auc_score(y, s)
hits = lambda s, k: int((rk(-s, method="min")[y] <= k).sum())

rule = {f: exn(f, "rand", 150) - exn(f, "ss", 150) for f in FOLD}
vote = {f: vfn(f, "ss", 150) for f in FOLD}
cons = np.mean([rk(v) for v in (*rule.values(), *vote.values())], 0)
print("== E19  Synergies found among the top-k ranked pairs (of 6; 153 pairs; random expectation = 6k/153)")
print(f"  {'score':34}" + "".join(f"{'k=' + str(k):>7}" for k in (5, 10, 15, 20, 35, 50)))
for name, s in (("paper rule, parameters only", rule["link-only"]), ("paper rule, topology + parameters", rule["topo-and-link"]),
                ("vote fraction, parameters only", vote["link-only"]), ("consensus of 6 scores", cons)):
    print(f"  {name:34}" + "".join(f"{hits(s, k):7}" for k in (5, 10, 15, 20, 35, 50)))
print(f"  {'random expectation':34}" + "".join(f"{6 * k / 153:7.1f}" for k in (5, 10, 15, 20, 35, 50)))

print("\n== E20  Fewer evolution runs: 50 runs (150 models) vs 150 runs (450 models)")
for n in (50, 150):
    comp = {}
    for f in FOLD:
        comp[f"vote {f}"] = vfn(f, "ss", n)
        try:
            comp[f"rule {f}"] = exn(f, "rand", n) - exn(f, "ss", n)
        except FileNotFoundError:
            pass
    c = np.mean([rk(v) for v in comp.values()], 0)
    single = ", ".join(f"{k.split()[0]} {k.split()[1][:4]} {auc(v):.2f}" for k, v in comp.items())
    print(f"  {n:3} runs: consensus of {len(comp)} scores AUC {auc(c):.3f} PR {average_precision_score(y, c):.3f} top-35 hits {hits(c, 35)}   | components: {single}")
