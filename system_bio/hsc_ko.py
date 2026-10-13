"""B: LOF/GOF scan of every node, HSC basin (synchronous, % of 2^21 states) vs paper claims."""
from hsc_model import *

# paper (Herrera 2024) claims, relative to WT HSC basin
PAPER = {("AKT", 1): "lost", ("FOXO3", 1): "up", ("Antiox", 1): "up", ("HIF1", 1): "up", ("mTOR", 1): "same",
         ("FOXO3", 0): "down", ("GATA2", 0): "lost", ("Meis1", 0): "lost", ("HIF1", 0): "lost",
         ("Cebpa", 0): "up", ("mTOR", 0): "same", ("OXPHOS", 0): "same"}


def scan(nxt):
    cycles, sizes = sync_basins(nxt)
    hsc = sum(z for c, z in zip(cycles, sizes) if len(c) == 1 and label(c[0]) == "HSC")
    return 100 * hsc / len(S), sum(len(c) == 1 for c in cycles), sum(len(c) > 1 for c in cycles)


wt, nf, nc = scan(NXT)
print(f"WT  HSC basin {wt:.4f}%  fixed {nf}  cyclic {nc}\n")
print(f"{'node':8}{'mut':5}{'HSC%':>9}{'fixed':>6}{'cyc':>5}  model   paper")
for g in names:
    for v in (0, 1):
        b, nf, nc = scan(clamp(NXT, g, v))
        rel = "lost" if b == 0 else "up" if b > wt + .01 else "down" if b < wt - .01 else "same"  # .01 %-point tolerance
        exp = PAPER.get((g, v))
        print(f"{g:8}{'GOF' if v else 'LOF':5}{b:9.4f}{nf:6}{nc:5}  {rel:6}  {exp or '':5}{'  OK' if exp == rel else '  <-- differs' if exp else ''}")
