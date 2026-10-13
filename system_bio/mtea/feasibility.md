# Benedicto 2025 feasibility check (2026-10-07)
Env: Python 3.11 venv tools/mtea (mteapy 0.x from GitHub, pydeseq2 0.5.4, cobra 0.32.1). Script: repro.py (de | tidee | tide; ALLGENES=1 = keep non-model genes in null).
Runtime (laptop, 1 core, ~1 GB free RAM): DESeq2 x5 ~11 min; TIDE-essential ~1-4 s/contrast; Human-GEM load ~16 s; TIDE ~8-14 s/contrast (1000 perms).

| cond | DEG mine | TIDE mine (model-genes null / all-genes null) | TIDE paper | TIDEe mine / ALL | TIDEe paper |
|---|---|---|---|---|---|
| TAKi (5Z) | 1782 | 14 / 20 | 25 | 10 / 10 | 14 |
| MEKi (PD) | 2400 | 29 / 17 | 17 | 12 / 12 | 15 |
| PI3Ki (PI103) | 1573 | 33 / 44 | 77 | 8 / 8 | 10 |
| PI3Ki-TAKi | 1986 | 42 / 26 | 55 | 7 / 5 | 11 |
| PI3Ki-MEKi | 3105 | 49 / 36 | 66 | 15 / 12 | 18 |

Paper: ~2000 DEG/condition. Task counts = p<0.025 (permutation). Mine: 143 tasks usable in TIDE-e (paper 187).
Spearman(score TIDE, score TIDE-e) per condition .53-.65 (paper pooled .93).
Deviations: pydeseq2 + apeglm shrink (not R DESeq2), gene symbols->Ensembl only for Human-GEM genes, 1000 perms (paper 10,000), design ~batch+condition (assumed).
