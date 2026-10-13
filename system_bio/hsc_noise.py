"""A: asynchronous stochastic simulation (MaBoSS-style, all rates 1) from the HSC state.
exp1 = paper's 'fluctuation': random initial state only.  exp2 = persistent noise: random spontaneous flips all along."""
from hsc_model import *

rng = np.random.default_rng(0)
FIX = S[NXT == S]
LAB = {int(s): label(s) for s in FIX}
HSC = next(int(s) for s in FIX if LAB[int(s)] == "HSC")
ACT = [g for g in names if HSC >> IDX[g] & 1]
INACT = [g for g in names if not HSC >> IDX[g] & 1 and g != "Oxygen"]  # Oxygen stays 0 (hypoxia)
BIT = 1 << np.arange(N)
KEEP = np.where(np.arange(N) == IDX["Oxygen"], 0.0, 1.0)  # Oxygen = environment: never noisy
T, n = 50.0, 20000  # paper: max_time=50


def run(s, eps_on=0.0, eps_off=0.0):
    """Gillespie: unstable node flips at rate 1; stable node flips spontaneously at eps_on (OFF->ON) / eps_off (ON->OFF)."""
    t = np.zeros(len(s))
    while True:
        on = (s[:, None] & BIT) > 0
        rate = ((NXT[s] ^ s)[:, None] & BIT > 0) + np.where(on, eps_off, eps_on) * KEEP
        R = rate.sum(1)
        t = t + rng.exponential(1 / np.maximum(R, 1e-12))
        go = (R > 0) & (t < T)
        if not go.any():
            return s
        k = (rate.cumsum(1) < (rng.random(len(s)) * R)[:, None]).sum(1)
        s = np.where(go, s ^ (1 << k), s)


def outcome(s):
    u, inv = np.unique(s, return_inverse=True)
    lab = np.array([LAB.get(int(x), "trans") for x in u])[inv]
    p = np.array([(lab == k).mean() for k in ("HSC", "MEP", "GMP", "CLP")])
    return p, -sum(x * np.log2(x) for x in p if x), (lab == "trans").mean()


def start(q_on, q_off):
    s = np.full(n, HSC)
    for g in INACT:
        s ^= (rng.random(n) < q_on).astype(np.int64) << IDX[g]
    for g in ACT:
        s ^= (rng.random(n) < q_off).astype(np.int64) << IDX[g]
    return s


def show(title, rows):
    print(f"\n{title}\n{'condition':22}{'HSC':>7}{'MEP':>7}{'GMP':>7}{'CLP':>7}{'trans':>7}{'H(bit)':>8}")
    for name, s in rows:
        p, h, tr = outcome(s)
        print(f"{name:22}" + "".join(f"{100 * x:7.1f}" for x in (*p, tr)) + f"{h:8.2f}")


conds = {"inactive nodes": (1, 0), "active nodes": (0, 1), "both": (1, 1)}
show("exp1: initial-state fluctuation only (paper Fig 6 design), q = flip prob per node",
     [(f"{c}, q={q}", run(start(q * a, q * b))) for c, (a, b) in conds.items() for q in (0.1, 0.3, 0.5)])
show("exp2: persistent noise from exact HSC state, eps = spontaneous flip rate per node (relaxation rate = 1)",
     [("none", run(np.full(n, HSC)))] +
     [(f"{c}, eps={e}", run(np.full(n, HSC), a * e, b * e))
      for c, (a, b) in {"OFF->ON flips": (1, 0), "ON->OFF flips": (0, 1), "both": (1, 1)}.items() for e in (0.005, 0.02, 0.05)])
