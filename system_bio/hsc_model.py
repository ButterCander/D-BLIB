"""Herrera et al. 2024 HSC Boolean model (21 nodes). State = int bitmask, bit i = names[i]."""
import re
import numpy as np

rules = {}
for l in open("hsc_rules.txt"):
    if "," in l and not l.startswith("targets"):
        n, r = l.split(",", 1)
        rules[n.strip()] = r.strip()
names = list(rules)
N = len(names)
IDX = {n: i for i, n in enumerate(names)}
S = np.arange(2 ** N, dtype=np.int64)

X = {n: (S >> i & 1).astype(bool) for i, n in enumerate(names)}
NXT = np.zeros_like(S)  # NXT[s] = synchronous successor of state s
for i, n in enumerate(names):
    NXT |= eval(re.sub(r"&\s*!", "& ~", rules[n]).replace("!", "~"), {}, X).astype(np.int64) << i

MARK = {"GATA2": "HSC", "GATA1": "MEP", "PU1": "GMP", "Gfi1": "CLP"}


def label(s):
    return "+".join(t for g, t in MARK.items() if s >> IDX[g] & 1) or "-"


def clamp(nxt, node, v):
    """LOF (v=0) / GOF (v=1): node fixed to v, like BoolNet fixGenes."""
    m, c = 1 << IDX[node], v << IDX[node]
    return nxt[(S & ~m) | c] & ~m | c


def sync_basins(nxt):
    """Synchronous attractors: ([cycle states], basin sizes)."""
    p = nxt
    for _ in range(N):  # p = nxt^(2^N): every state lands on its attractor
        p = p[p]
    cid, cycles = np.full(len(S), -1), []
    for r in np.unique(p):
        if cid[r] < 0:
            c, x = [int(r)], nxt[r]
            while x != r:
                c.append(int(x))
                x = nxt[x]
            cid[c] = len(cycles)
            cycles.append(c)
    return cycles, np.bincount(cid[p], minlength=len(cycles))
