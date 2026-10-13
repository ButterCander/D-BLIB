import sys, json
sys.argv = ["x"]
src = open("eflux_grid.py").read()
head = src.split("MEAS = {")[0]
ns = {}; exec(compile(head, "e", "exec"), ns)
import numpy as np
from scipy import stats
E = ns["E"]; solve = ns["solve"]; CONDS = ns["CONDS"]
MEAS = {"PI103": 62.5, "PD0325901": 47.0, "5Z": 59.0, "PI103+5Z": 92.0, "PI103+PD": 78.0}
out = json.load(open("out/eflux_grid.json"))
for fl in (0.05, 0.2):
    for vmax in (0.3, 1.0, 3.0):
        g = np.array([solve(c, vmax, fl) for c in E.columns]).reshape(4, 6); r = g / g[:, [0]]
        pred = {c: 100 * (1 - np.nanmean(r[:, i])) for i, c in enumerate(CONDS) if i > 0}
        a = [pred[c] for c in MEAS]; b = [MEAS[c] for c in MEAS]
        out.append(dict(vmax=vmax, floor=fl, growth=g.tolist(), pred=pred, rho=float(stats.spearmanr(a, b)[0]), conc=0))
        print(vmax, fl, {k: round(v, 1) for k, v in pred.items()}, flush=True)
json.dump(out, open("out/eflux_grid.json", "w")); print("done")
