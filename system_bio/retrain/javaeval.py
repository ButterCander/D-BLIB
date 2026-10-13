"""Evaluate Java (Gitsbe + Drabme) runs: paper rule = random-minus-calibrated ensemble-wise Bliss excess, and calibrated vote fraction.
usage: python javaeval.py <label> <ss_run_name> <rand_run_name>   (folders under retrain/java_runs)"""
import glob
import os
import sys
import warnings
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
pk = lambda s: frozenset(s.replace("[", "").replace("]", "").split("-"))
OBS = {frozenset(l.strip().split("~")) for l in open(os.path.join(HERE, "..", "ags", "data", "observed_synergies_cascade_2.0")) if l.strip()}


def find(run, kind):
    hits = glob.glob(os.path.join(HERE, "java_runs", run, "**", f"*_{kind}_synergies.tab"), recursive=True)
    return pd.read_csv(hits[0], sep="\t", index_col=0) if hits else None


label, ss, rd = sys.argv[1:4]
E = {k: find(r, "ensemblewise") for k, r in (("ss", ss), ("rand", rd))}
Mw = {k: find(r, "modelwise") for k, r in (("ss", ss), ("rand", rd))}
if any(v is None for v in (*E.values(), *Mw.values())):
    sys.exit(f"missing output tables: ensemblewise {[k for k, v in E.items() if v is None]}, modelwise {[k for k, v in Mw.items() if v is None]}")
for d in (*E.values(), *Mw.values()):
    d.index = [pk(i) for i in d.index]
pairs = list(Mw["ss"].index)
y = np.array([p in OBS for p in pairs])
ex = lambda k: E[k].iloc[:, 0].reindex(pairs).fillna(0).values
vote = (Mw["ss"].Synergies / (Mw["ss"].Synergies + Mw["ss"]["Non-synergies"])).fillna(0).values
rule = ex("rand") - ex("ss")
n_models = int((Mw["ss"].Synergies + Mw["ss"]["Non-synergies"]).max())
print(f"{label}: {len(pairs)} pairs, {int(y.sum())} positives, up to {n_models} model evaluations per pair")
print(f"  paper rule AUC {roc_auc_score(y, rule):.3f} (PR {average_precision_score(y, rule):.3f}) | vote fraction AUC {roc_auc_score(y, vote):.3f} (PR {average_precision_score(y, vote):.3f})")
