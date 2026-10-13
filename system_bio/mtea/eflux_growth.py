"""E-Flux (expression-scaled capacity) FBA on Human-GEM with a Ham's-F12-like medium; GSE285616 24 h samples.
Question: does the predicted biomass flux of each treatment (relative to the batch-matched control) track the measured 48 h growth inhibition
(Flobak 2015: PI103+5Z 92%, PI103+PD 78%, PI103 52-73%, 5Z 59%, PD 47%)?   usage: python eflux_growth.py [Vmax ...]"""
import ast, sys, time, json
import numpy as np, pandas as pd
from cobra.io import read_sbml_model
D = "mteapy/src/mteapy/data"
CONDS = ["control", "PI103", "PD0325901", "5Z", "PI103+5Z", "PI103+PD"]
VMAX = [float(v) for v in sys.argv[1:]] or [10.0]
model = read_sbml_model(f"{D}/HumanGEM.xml.gz")
BIO = "MAR13082"
gm = pd.read_csv(f"{D}/HumanGEM_genes.tsv", sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]

# ---- expression: size-factor normalised counts, gene IDs of Human-GEM
cnt = pd.read_csv("GSE285616_counts.csv.gz", index_col="Symbol"); cnt = cnt[~cnt.index.duplicated()]
x = cnt.values.astype(float); pos = (x > 0).all(axis=1)
sf = np.median(x[pos] / np.exp(np.log(x[pos]).mean(axis=1, keepdims=True)), axis=0)
E = pd.DataFrame(x / sf, index=cnt.index, columns=cnt.columns)
E = E[E.index.isin(gm.index)]; E.index = gm[E.index].values
samples = [(c, b) for b in range(4) for c in range(6)]    # column order: batch-major, 6 conditions per batch
FLOOR = 1.0

# ---- GPR evaluation (AND = min, OR = sum)
def ev(node, e):
    if isinstance(node, ast.Name): return e.get(node.id, FLOOR)
    if isinstance(node, ast.BoolOp):
        v = [ev(n, e) for n in node.values]
        return min(v) if isinstance(node.op, ast.And) else sum(v)
    if isinstance(node, ast.Expression): return ev(node.body, e)
    raise ValueError(type(node))
gpr = {r.id: r.gpr.body for r in model.reactions if r.gpr is not None and r.gpr.body is not None}
rxn_ids = list(gpr)
print(f"reactions with GPR: {len(gpr)} of {len(model.reactions)}", flush=True)
def rxn_expr(col):
    e = E[col].to_dict()
    return pd.Series({rid: ev(node, e) for rid, node in gpr.items()})
RE = pd.DataFrame({col: rxn_expr(col) for col in E.columns})
ref = RE[[c for c in E.columns if (list(E.columns).index(c) % 6) == 0]].mean(axis=1).median()
print(f"reference expression (median of control reactions): {ref:.1f}", flush=True)

# ---- medium
AA = "alanine arginine asparagine aspartate cysteine glutamate glutamine glycine histidine isoleucine leucine lysine methionine phenylalanine proline serine threonine tryptophan tyrosine valine".split()
ex_by_name = {}
for r in model.exchanges:
    mets = list(r.metabolites)
    if len(mets) == 1: ex_by_name[mets[0].name.lower()] = r.id
MED = {"glucose": 10, "pyruvate": 1, "glutamine": 3, "choline": .5, "inositol": .5, "biotin": .05, "folate": .05, "nicotinamide": .05, "pantothenate": .05,
       "pyridoxine": .05, "riboflavin": .05, "thiamin": .05, "lipoic acid": .05, "hypoxanthine": .1, "thymidine": .1, "putrescine": .1, "linoleate": .1,
       "cholesterol": .1, "palmitate": .1, "oleate": .1, "ethanolamine": .1, "o2": 30}
for a in AA: MED.setdefault(a, 1.0)
OPEN = ["h2o", "co2", "h+", "na+", "k+", "ca2+", "mg2+", "fe2+", "fe3+", "cu2+", "zinc", "sulfate", "pi", "chloride", "hco3-", "nh3"]
missing = [n for n in list(MED) + OPEN if n not in ex_by_name]
print("medium components without exchange:", missing, flush=True)
default = {}
for r in model.reactions: default[r.id] = (r.lower_bound, r.upper_bound)
exch = {r.id for r in model.exchanges}
def apply_medium():
    for rid in exch:
        model.reactions.get_by_id(rid).bounds = (0.0, 1000.0)
    for n, u in MED.items():
        if n in ex_by_name: model.reactions.get_by_id(ex_by_name[n]).bounds = (-u, 1000.0)
    for n in OPEN:
        if n in ex_by_name: model.reactions.get_by_id(ex_by_name[n]).bounds = (-1000.0, 1000.0)
apply_medium()
from cobra import Reaction
_b = model.reactions.get_by_id(BIO)
_x = Reaction("BIOX"); _x.name = "biomass without lipid and cofactor/vitamin pools"; _x.bounds = (0.0, 1000.0)
model.add_reactions([_x])
_x.add_metabolites({m: c for m, c in _b.metabolites.items() if m.id not in ("MAM10012c", "MAM10014c")})
default["BIOX"] = (0.0, 1000.0)
model.reactions.get_by_id(BIO).bounds = (0.0, 0.0)
BIO = "BIOX"
model.objective = BIO

def solve(col, vmax, shuffle=None):
    re = RE[col] if shuffle is None else pd.Series(RE[col].values[shuffle], index=RE.index)
    cap = vmax * re / ref
    for rid in rxn_ids:
        if rid in exch or rid == BIO: continue
        lb, ub = default[rid]; c = float(cap[rid])
        nlb = max(lb, -c) if lb < 0 else lb
        nub = min(ub, c) if ub > 0 else ub
        if nlb > nub: nlb = nub
        model.reactions.get_by_id(rid).bounds = (nlb, nub)
    s = model.optimize()
    return s.objective_value if s.status == "optimal" else np.nan

res = []
for vmax in VMAX:
    t0 = time.time()
    g = []
    for i, col in enumerate(E.columns):
        g.append(solve(col, vmax))
    g = np.array(g).reshape(4, 6)       # batch x condition
    ratio = g / g[:, [0]]
    print(f"\nVmax={vmax:g}  ({time.time()-t0:.0f}s)  growth per condition (mean over batches): " + ", ".join(f"{c} {g[:, k].mean():.3f}" for k, c in enumerate(CONDS)), flush=True)
    print("   ratio to batch-matched control (mean +- sd): " + ", ".join(f"{c} {ratio[:, k].mean():.3f}+-{ratio[:, k].std(ddof=1):.3f}" for k, c in enumerate(CONDS[1:], 1)), flush=True)
    res.append(dict(vmax=vmax, growth=g.tolist(), ratio=ratio.tolist()))
json.dump(res, open("out/eflux_results.json", "w"))
