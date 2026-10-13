"""Link-only: batch-to-batch variation and the pooled larger ensembles (all GA seeds distinct)."""
import os
import pickle
import warnings
import numpy as np
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
rk = rankdata
OBS = [frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()]
ld = lambda name: pickle.load(open(os.path.join(HERE, "out", f"balanced_{name}.pkl"), "rb"))


def part(name, skip=0):
    d = ld(name)
    return {"responses": d["responses"][skip:], "fitness": d["fitness"][skip:], "min_output": d["min_output"], "max_output": d["max_output"]}


def join(parts):
    out = dict(parts[0])
    out["responses"] = [x for p in parts for x in p["responses"]]
    out["fitness"] = [x for p in parts for x in p["fitness"]]
    return out


def scores(d):
    lo, hi = d["min_output"], d["max_output"]
    n = lambda v: (v - lo) / (hi - lo)
    resp = d["responses"]
    pairs = sorted({k for r in resp for k in r if "-" in k})
    ens, vote = [], []
    for p in pairs:
        a, b = p.split("-")
        m = lambda k: np.mean([r[k] for r in resp if r.get(k) is not None] or [np.nan])
        ens.append(n(m(p)) - n(m(a)) * n(m(b)))
        s = [n(r[p]) - n(r[a]) * n(r[b]) < 0 for r in resp if all(r.get(k) is not None for k in (a, b, p))]
        vote.append(np.mean(s) if s else 0.0)
    return pairs, np.nan_to_num(np.array(ens)), np.array(vote)


CH = (100, 112, 124)
SS = {"seeds 12-23": part("ss_12"), "seeds 24-47": part("ss_24"), "seeds 48-71": part("ss_36", 36), **{f"seeds {b}-{b + 11}": part(f"ss_12_b{b}") for b in CH}}
RD = {"seeds 12-23": part("rand_12"), "seeds 24-47": part("rand_24"), **{f"seeds {b}-{b + 11}": part(f"rand_12_b{b}") for b in CH}}
pairs, e_r_all, v_r_all = scores(join(list(RD.values())))
y = np.array([frozenset(p.split("-")) in OBS for p in pairs])
auc = lambda s: roc_auc_score(y, s)
print("== calibrated ensembles by seed batch, rule computed against the pooled random-model ensemble")
rows = []
for name, d in SS.items():
    _, e, v = scores(d)
    rule = e_r_all - e
    rows.append((auc(rule), auc(v), auc(rk(rule) + rk(v))))
    print(f"  {name:14} models {len(d['fitness']):3} | rule {rows[-1][0]:.3f} | vote {rows[-1][1]:.3f} | 2-score consensus {rows[-1][2]:.3f}")
r = np.array(rows)
print(f"  spread between batches (sd): rule {r[:, 0].std():.3f}, vote {r[:, 1].std():.3f}, consensus {r[:, 2].std():.3f}")
full = join(list(SS.values()))
_, e, v = scores(full)
rule = e_r_all - e
c = rk(rule) + rk(v)
print(f"\n== POOLED link-only: {len(full['fitness'])} calibrated models, {len(join(list(RD.values()))['fitness'])} random models")
print(f"  rule {auc(rule):.3f} | vote {auc(v):.3f} | 2-score consensus {auc(c):.3f} (PR {average_precision_score(y, c):.3f}, top-35 {int((rk(-c, method='min')[y] <= 35).sum())}/6)")
print("  published tables (450 models): rule 0.686 | vote 0.831 | 2-score consensus 0.827")
