"""Re-train CASCADE 2.0 Boolean model ensembles with pyDrugLogics and simulate every drug / drug pair, per model.
usage: python pipeline.py <balanced|topology|mixed> <ss|rand> <runs> [cores]   -> retrain/out/<mt>_<kind>_<runs>.pkl"""
import logging
import os
import pickle
import sys
import time
import warnings
from joblib import Parallel, delayed
from pydruglogics.input.ModelOutputs import ModelOutputs
from pydruglogics.input.Perturbations import Perturbation
from pydruglogics.input.TrainingData import TrainingData
from pydruglogics.model.BooleanModel import BooleanModel
from pydruglogics.model.Evolution import Evolution
from pydruglogics.model.InteractionModel import InteractionModel

warnings.filterwarnings("ignore")
logging.disable(logging.CRITICAL)
BACKEND = "threading" if os.environ.get("THREADS") else "loky"   # threads share memory (a loky worker costs ~150 MB); clingo calls release the GIL
if BACKEND == "threading":
    import joblib
    import pydruglogics.model.Evolution as _E
    _E.Parallel = lambda *a, **k: joblib.Parallel(*a, **{**k, "backend": "threading"})
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "ags_cascade_2.0") + os.sep


class SerialEvolution(Evolution):
    """pyDrugLogics nests joblib pools (runs x population) -> hundreds of processes and out-of-memory; evaluate the population serially instead."""

    def calculate_fitness(self, ga_instance, solutions, solution_idx):
        return [self._calculate_fitness_for_solution(s) * 100 for s in solutions]


def read(f):
    return [l.rstrip("\n") for l in open(D + f) if l.strip() and not l.startswith("#")]


def simulate(model, perts, mo):
    out = {}
    for p in perts:
        m = model.clone()
        m.add_perturbations(p)
        m.calculate_attractors(m.attractor_tool, m.attractor_type)
        out["-".join(d["name"] for d in p)] = m.calculate_global_output(mo, False)
    return out


def main(mt, kind, runs, cores):
    mo = ModelOutputs(input_dictionary={"RSK_f": 1.0, "CCND1": 1.0, "MYC": 1.0, "FOXO_f": -1.0, "CASP8": -1.0, "CASP9": -1.0})
    drug_data = [[l.split("\t")[0], l.split("\t")[2:], "inhibits"] for l in read("drugpanel")]
    perts = Perturbation(drug_data=drug_data, perturbation_data=[l.split("\t") for l in read("perturbations")])
    im = InteractionModel(interactions_file=D + "network.sif", model_name="cascade2", remove_inputs=True)
    bm = BooleanModel(model=im, mutation_type=mt, attractor_tool="mpbn", attractor_type="stable_states")
    if kind == "ss":
        td = TrainingData(input_file=D + "training")
    elif kind == "aug":  # new: 24 literature nodes + expression-derived extra nodes (retrain/data/extra_annotations.tsv)
        lines = read("training")
        nodes = lines[[l.lower() for l in lines].index("response") + 1].split("\t")
        extra = [f"{r.split(chr(9))[0]}:{r.split(chr(9))[1]}" for r in open(os.path.join(HERE, "data", "extra_annotations.tsv")).read().splitlines()[1:] if r.strip()]
        td = TrainingData(observations=[(nodes + extra, 1.0)])
    elif kind.startswith("ctl"):  # negative control: 24 literature nodes + 8 RANDOM other nodes with RANDOM 0/1 states (ctl<seed>)
        import pickle as _pk
        import random
        lines = read("training")
        nodes = lines[[l.lower() for l in lines].index("response") + 1].split("\t")
        nd = _pk.load(open(os.path.join(HERE, "data", "nodes.pkl"), "rb"))
        rng_ = random.Random(int(kind[3:]))
        pick = rng_.sample([n for n in nd["nodes"] if n not in nd["trained"]], 8)
        td = TrainingData(observations=[(nodes + [f"{n}:{rng_.randint(0, 1)}" for n in pick], 1.0)])
    elif kind.startswith("sub"):  # sub<k>_<seed>: random subset of k of the 24 steady-state nodes (calibration-data ablation)
        import random
        k, sd = kind[3:].split("_")
        lines = read("training")
        nodes = lines[[l.lower() for l in lines].index("response") + 1].split("\t")
        td = TrainingData(observations=[(random.Random(int(sd)).sample(nodes, int(k)), 1.0)])
    elif kind == "hyb":  # new: steady-state match (24 nodes) AND proliferation, equal weights, in one fitness
        lines = read("training")
        nodes = lines[[l.lower() for l in lines].index("response") + 1].split("\t")
        td = TrainingData(observations=[(nodes, 1.0), (["globaloutput:1"], 1.0)])
    else:
        td = TrainingData(observations=[(["globaloutput:1"], 1.0)])
    t0 = time.time()
    evo = SerialEvolution(boolean_model=bm, model_outputs=mo, training_data=td,
                    ga_args={"num_generations": 20, "num_parents_mating": 3, "mutation_num_genes": 3, "fitness_batch_size": 20},
                    ev_args={"num_best_solutions": 3, "num_of_runs": runs, "num_of_cores": cores, "num_of_init_mutation": 10})
    base = int(os.environ.get("SEED_BASE", runs))   # GA seeds are base..base+runs-1; batches with different (runs) can overlap, so allow an explicit base
    evo._global_seed = base
    import numpy as _np
    _np.random.seed(base)
    models = evo.run()
    t1 = time.time()
    res = Parallel(n_jobs=cores, backend=BACKEND)(delayed(simulate)(m, perts.perturbations, mo) for m in models)
    t2 = time.time()
    out = {"mt": mt, "kind": kind, "runs": runs, "fitness": [m.fitness / 100 for m in models], "names": [m.model_name for m in models],
           "responses": res, "min_output": mo.min_output, "max_output": mo.max_output, "train_s": t1 - t0, "sim_s": t2 - t1}
    os.makedirs(os.path.join(HERE, "out"), exist_ok=True)
    out["seed_base"] = base
    suffix = f"_b{base}" if "SEED_BASE" in os.environ else ""
    pickle.dump(out, open(os.path.join(HERE, "out", f"{mt}_{kind}_{runs}{suffix}.pkl"), "wb"))
    print(f"{mt} {kind} runs={runs}: {len(models)} models, train {t1 - t0:.0f}s, simulate {t2 - t1:.0f}s, fitness mean {sum(out['fitness']) / len(models):.3f} max {max(out['fitness']):.3f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]) if len(sys.argv) > 4 else 3)
