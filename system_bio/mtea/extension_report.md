# Extension analyses on Benedicto 2025 (not in the slides)

Data: GEO GSE285616 counts; authors' Supplementary Data 3, 4, 6 (Nature, static-content.springer.com, 0.2 MB, saved in `si/`).
Code: `repro.py` (DE, TIDE, TIDE-essential; `NPERM`, `ALLGENES` env), `inter.py` (interaction contrasts), `camera_like.py`, `taskgenes.py`. Outputs in `out/`.
Convention: two-sided p = 2 × the paper's tail p; BH = Benjamini–Hochberg across all tasks of one condition.

## 1. BH on my reproduction (pydeseq2 + MTEApy, 10,000 permutations, null = all genes)
| condition | TIDE p<0.025 → BH q<0.05 | TIDE-essential |
|---|---|---|
| TAKi | 22 → 0 | 10 → 0 |
| MEKi | 16 → 0 | 13 → 0 |
| PI3Ki | 44 → 2 | 8 → 0 |
| PI3Ki+TAKi | 29 → 3 | 5 → 0 |
| PI3Ki+MEKi | 34 → 3 | 16 → 0 |
Survivors: serine synthesis, glycine synthesis (PI3Ki-containing conditions), UDP-GalNAc synthesis, cysteine degradation.
Asp → Asn (ASNS) does not survive: a single-gene task has p of about rank / #genes, so it cannot get small enough for 187 tests.
With 1,000 permutations the p resolution (0.001) limits BH when few tasks reach the minimum p; 10,000 permutations give finer resolution (a precision gain, not a bias correction). The table above uses the all-genes null; the model-gene null (MTEApy default) gives different counts, see section 11.

## 2. BH on the authors' own p-values (Supplementary Data 3)
| condition | p<0.025 | BH q<0.05 | q<0.10 |
|---|---|---|---|
| TAKi | 31 | 0 | 0 |
| MEKi | 19 | 0 | 0 |
| PI3Ki | 78 | 59 | 69 |
| PI3Ki+TAKi | 56 | 3 | 32 |
| PI3Ki+MEKi | 64 | 27 | 43 |
NB: these BH counts test the combination against CONTROL (is the task changed at all?), not synergy (combination against the single drugs).
Authors' synergy calls (Supplementary Data 4), BH on the combination-vs-control p-values:
- PI3Ki+MEKi: 48 calls → 23 at q<0.05, 36 at q<0.10.
- PI3Ki+TAKi: 14 calls → 1 at q<0.05, 8 at q<0.10.
- Strong synergy (4 + 4 calls): every call has a single-drug p between 0.025 and 0.10 (a near-miss of the threshold); combination q = 0.127–0.158 (spermidine 0.127, CoQ 0.135, dTTP 0.135, aspartate-from-glutamine 0.135, thioredoxin 0.158); none passes q<0.10.
- Weak synergy: median margin over the best single drug 0.076–0.081 log2FC; 60–68% of calls below 0.10.
Supplement inconsistencies: sheet counts 31/19/78/56/64 vs text 30/17/77/55/66; sheets "TAKi" and "MEKi" carry swapped header cells; 48 synergy rows vs "49" in the text.
My scores vs authors' scores (187 tasks): r = 0.95 (PI3Ki), 0.94 (PI3Ki+TAKi), 0.78 (PI3Ki+MEKi), 0.82 (TAKi), 0.40 (MEKi; Spearman 0.75). My p-values are more conservative (median p 0.11–0.20 vs 0.05–0.13).

## 3. Interaction model (same counts, `~ batch + condition`, custom contrasts, unshrunk log2FC, 3,000 permutations, null = all genes)
Dose-scaled expectation: AB − 0.5·A − 0.5·B (combination is half dose of each drug; assumes linear dose–response).
| combination | method | p<0.025 | BH q<0.05 | overlap with authors' synergy calls (raw / BH) |
|---|---|---|---|---|
| PI3Ki+MEKi | TIDE | 76 (75 negative) | 41 | 45 / 32 of 48 |
| PI3Ki+MEKi | TIDE-essential | 24 | 0 | 10 / 0 |
| PI3Ki+TAKi | TIDE | 27 (20 negative) | 0 | 7 / 0 of 14 |
| PI3Ki+TAKi | TIDE-essential | 9 | 0 | 1 / 0 |
Full-dose additive contrast (AB − A − B) finds almost nothing (expected for a half-dose combination).
Strongest PI3Ki+MEKi interaction tasks: serine synthesis (−0.81), glycine synthesis (−0.75), then amino-acid degradation tasks (−0.13 to −0.23).
Authors' strong-synergy examples under the interaction test: dTTP q = 0.037, spermidine q = 0.071, CoQ q = 0.126, aspartate-from-glutamine q = 0.129.
Reading: PI3Ki+MEKi non-additivity is supported (broad, mostly amino-acid related); PI3Ki+TAKi is not. A saturating dose–response would make the half-dose expectation too low and inflate this result.

## 4. CAMERA-like, correlation-aware competitive test (my simplified implementation, Wald z of DESeq2, residual correlations from 15 residual df)
Median residual correlation inside a set: 0.18 (essential-gene sets, median size 6) and 0.10 (task gene sets, median size 17) → variance inflation 1.6 and 2.6.
Task gene sets, two-sided p<0.05:
| condition | independence | limma default rho = 0.01 | estimated rho | BH q<0.05 (estimated rho) |
|---|---|---|---|---|
| TAKi | 63 | 59 | 12 | 0 |
| MEKi | 49 | 37 | 15 | 0 |
| PI3Ki | 84 | 61 | 11 | 0 |
| PI3Ki+TAKi | 44 | 34 | 13 | 2 |
| PI3Ki+MEKi | 76 | 56 | 12 | 0 |
Essential-gene sets: gene-shuffle-significant tasks (4–13 per condition) mostly lose significance under the estimated correlation (2–9 lost per condition).

## Caveats
- My reproduction is more conservative than the authors' own p-values; conclusions in section 2 use the authors' numbers.
- Interaction contrasts are not shrunk (pydeseq2 shrinkage needs a named coefficient).
- Section 4 is a simplified CAMERA, not limma; sets differ from TIDE's GPR-projected reaction scores; correlations from 24 samples are noisy.
- Everything here is computational; no experimental validation.

## 5. Highest-single-agent test (combination vs each single drug; same counts, unshrunk contrasts, 3,000 permutations, null = all genes)
A task counts as better than the best single drug only if the combination is significantly more extreme than BOTH singles (intersection-union p = max of the two, then BH).
| combination | method | vs PI3Ki only (p<.025) | vs partner only (p<.025) | both, same sign | BH q<.05 | authors' synergy calls supported (raw / BH) |
|---|---|---|---|---|---|---|
| PI3Ki+MEKi | TIDE | 14 | 92 | 2 | 0 | 2 / 0 of 48 |
| PI3Ki+MEKi | TIDE-essential | 8 | 33 | 0 | 0 | 0 / 0 |
| PI3Ki+TAKi | TIDE | 16 | 69 | 3 | 0 | 3 / 0 of 14 |
| PI3Ki+TAKi | TIDE-essential | 4 | 22 | 0 | 0 | 0 / 0 |
Examples (PI3Ki+MEKi, score of combination minus single): serine synthesis vs PI3Ki +0.11 (p = 0.44), vs MEKi −1.02; Asp → Asn vs PI3Ki 0.00, vs MEKi −3.35; spermidine synthesis vs PI3Ki −0.43 (p = 0.052), vs MEKi −0.37 (p = 0.034).
PI3Ki+TAKi, Asp → Asn: combination is significantly LESS extreme than PI3Ki alone (+1.48, p = 0.012).
Reading: the half-dose combination looks like full-dose PI3Ki alone, not better. The authors' weak-synergy calls (margin about 0.08) are not distinguishable from "equal to PI3Ki".
The design has no half-dose single-drug arms, so "equal effect at half dose" cannot be separated from a saturating PI3Ki dose–response.

## Combined reading
| evidence | PI3Ki+MEKi | PI3Ki+TAKi |
|---|---|---|
| authors' calls with combination-vs-control p significant after BH (NOT a synergy test) | 23 of 48 (q<.05) | 1 of 14 |
| interaction vs half-dose expectation | 41 tasks (BH), 32 of 48 calls | 0 |
| better than best single drug (HSA) | 0 | 0 |
| strong-synergy examples | q .13–.16, near-miss thresholds | same |
Conclusion: PI3K inhibition is the robust effect; the combination matches it at half dose; metabolic "synergy" is not established by this design.

## 6. Verification of the analyses above (`validate_checks.py`, `validate_interaction.py`, `nullde.py`)
| check | result | consequence |
|---|---|---|
| contrast algebra: AB−A, AB−B, AB−.5A−.5B, AB−A−B equal the combination of single-condition log2FC | max difference 1e-15 | contrasts are right |
| BH: scipy vs manual step-up | max difference 1e-16 | implementation is right |
| p convention: two-sided 2p with q<.05 vs tail p with q<.05 (authors, PI3Ki+MEKi) | 27 vs 43 tasks | two-sided is the correct match to the paper's p<0.025 |
| full pipeline on shuffled labels (pydeseq2, condition shuffled within batch, 5 contrasts) | 0 genes with padj<.05, 0 tasks | DE step and masking give no false positives; but the null is degenerate (nothing to score) |
| set-level tests on shuffled labels (150 null contrasts, unmasked Wald z, real gene sets) | gene-label permutation (MTEApy logic): false-positive rate 0.39 (task sets) / 0.22 (essential sets) at p<.05, about 35 / 8 BH rejections per null contrast; CAMERA-like with estimated rho: 0.037 / 0.026, 0.02 / 0.00 rejections | gene-shuffle p-values are anti-conservative for correlated sets; CAMERA-like is calibrated |
| spike-in power (known shift added to one task's genes, 60 trials) | CAMERA-like detects 0.1 / 0.2 / 0.3 / 0.5 log2 shifts at p<.05 in 23 / 57 / 80 / 95 % (BH: 3 / 25 / 42 / 77 %); gene-shuffle 63 / 80 / 97 / 98 % | CAMERA-like is conservative but not blind to shifts of 0.3 or more |
| detectable task score for the HSA contrast (set sizes 6 / 17 / 50) | |score| ≥ 0.60 / 0.35 / 0.22 at p<.05 | authors' margins (about 0.08) are far below what this design can resolve |
| calibrated re-test of the contrasts (CAMERA-like, 120 task sets) | dose-scaled interaction PI3Ki+MEKi: 40 tasks at p<.05 (38 negative), 0 at BH q<.05 (smallest q = 0.075, serine and glycine synthesis); PI3Ki+TAKi: 11 / 0; AB vs PI3Ki alone: 8 tasks, all positive (combination less extreme), 0 negative; AB vs MEKi alone: 93, BH 87 | interaction support is suggestive, not established |

Corrections to earlier statements:
- "The half-dose combination equals PI3Ki alone" → "no difference is detectable"; the design cannot resolve differences below about 0.2–0.6 log2 per task, and the calibrated test even shows the combination slightly LESS extreme than PI3Ki alone for 8 tasks.
- Interaction result (41 tasks) used gene-shuffle p-values; under the calibrated test no task passes FDR 5%.
- Counts based on the authors' gene-shuffle p-values (section 2) are upper bounds for the same reason.

## 7. Checking against measured drug efficacy (Flobak et al. 2015, PLoS Comput Biol, S1 Text; AGS real-time growth assays)
Source file: PLOS supplementary S1 Text (35 pages, 3.2 MB) saved as `si/flobak_S1.bin`; text in `si/flobak_S1.txt`; growth-curve images extracted to `si/fl_p18_img*.png`.
Design in that study: single drugs at GI50 and two 2-fold dilutions; combinations at 0.5x GI50 of each drug (same concentrations as the RNA-seq study); growth 48 h after adding drugs; synergy declared when the half-dose combination beats each single drug at full GI50; Chou–Talalay CI reported "well below 0.5" for PI3K–MEK, TAK1–PI3K, MEK–AKT, TAK1–AKT.
Digitized end of assay (cell index 70 h; drug added at 21 h, index normalized to 1; one representative replicate, about ±0.1):
| experiment | control | single A | single B | combination (half doses) |
|---|---|---|---|---|
| PI103 + PD0325901 | +3.25 | PI103 +1.55 (52% inhibited) | PD +1.71 (47%) | +0.73 (78%) |
| PI103 + (5Z)-7-oxozeaenol | +3.50 | PI103 +0.96 (73%) | 5Z +1.45 (59%) | +0.27 (92%) |
Growth-level HSA is clearly met (combination increment 0.47x and 0.28x of the best single drug). Both combinations also exceed the Bliss expectation computed from the half-dose single effects.
Dose–response shape (fitted curves, S1 Fig 3, 7, 9, read by eye): half-GI50 gives about 0.7–0.8 of the full-GI50 growth inhibition (PI103 ≈ 35% vs 50%; PD0325901 ≈ 29% vs 36%) — concave, not linear.
Interaction test with weight w for the expected half-dose effect (expected AB = w·A + w·B; CAMERA-like, 120 task sets):
| combination | w = 0.5 | 0.6 | 0.7 | 0.7/0.8 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|---|
| PI3Ki+MEKi tasks p<.05 (BH q<.05) | 34 (0) | 38 (0) | 33 (1) | 31 (1) | 15 (1) | 5 (0) | 7 (0) |
| PI3Ki+TAKi tasks p<.05 (BH q<.05) | 13 (0) | 3 (1) | 2 (1) | 2 (1) | 2 (1) | 3 (0) | 5 (0) |
With the growth-informed weights the PI3Ki+MEKi signal shrinks and is not established; PI3Ki+TAKi has no interaction signal.
Proliferation check on the 24 h RNA-seq (Tirosh S + G2/M genes, 94 genes; mean log2FC vs control): PI3Ki +0.31, TAKi +0.11, PI3Ki+TAKi +0.25, MEKi −0.42, PI3Ki+MEKi −0.40; PI3Ki+MEKi vs PI3Ki alone −0.70 but vs MEKi alone +0.02.
The cell-cycle signature does not follow the 48 h growth ranking (growth inhibition: PI3Ki+TAKi 92% > PI3Ki+MEKi 78% > PI3Ki 52–73% ≈ TAKi 59% > MEKi 47%); the strongest-inhibited PI3Ki+TAKi condition has cell-cycle genes above control.
Reading: the drugs and combinations have real, measured growth synergy; the 24 h transcriptome (and the metabolic task scores derived from it) does not track it, so "metabolic synergy" labels cannot be validated or used as a proxy for efficacy.
Caveats: digitized from published figures (single replicate), different experiment and time point (24 h RNA vs 48 h growth), CI values not tabulated numerically in the supplement.

## 8. Independent, dose-matched, time-resolved RNA-seq (Tsirvouli et al. 2024 repository; `tsirvouli_analysis.py`)
Source: github.com/Eirinits/AGS_synergies, Transcriptomics/Data (`gene_counts_flt_norm_wCTRL.csv` 24 MB, `gene_info.tsv` 3.7 MB, `sample_info_to_use.tsv`; saved in `tsirvouli/`). 94 samples = 6 treatments (DMSO, 5Z, PD, PI, PI+5Z, PI+PD) × 1, 2, 4, 8, 24 h × 3 replicates (+ untreated 0 h), 13,565 genes, log2-normalised values; analysed per time point with a cell-means linear model and limma-style variance moderation.
Dosing: the Methods of the preprint give PI103 0.70 µM, 5Z 0.50 µM, PD0325901 0.035 µM for single drugs and combinations (full doses). Scale: `gene_counts_flt_norm_wCTRL.csv` is VSN (glog2) output (official script `tr_limma_per_Timepoint.R`: `vsn::predict`), so the 'log2FC' values here are VSN-scale differences, compressed toward 0 for low-expressed genes (matrix minimum 5.3; 23.5% of genes have mean < 7); official limma-voom logFC was NOT used (see section 11).
| check | result |
|---|---|
| dose diagnostic: slope of LFC(combination) on LFC(single) for genes changed by one drug only | this dataset 0.76–0.95 at 1–8 h (0.49 for PI at 24 h); GSE285616 at 24 h 0.47–0.75 → combinations in GSE285616 carry about half the single-drug effect, consistent with the half-dose design (selection on the noisy x variable biases all slopes downward) |
| 24 h replication of gene-level log2FC vs DMSO (11,171 shared genes) | PI r = 0.60, PD 0.79, 5Z 0.54, PI+PD 0.77, PI+5Z 0.72 |
| key genes at 24 h (log2FC): ASNS | PI −1.8, PD −1.1, 5Z −0.2, PI+PD −3.8, PI+5Z −2.7 (GSE285616: PI −5.7, PD −1.9, PI+PD −5.4): same direction, smaller PI effect here |
| task-set tests (109 sets, CAMERA-like, estimated rho) | at every time point 0–9 sets at p<.05 and 0–3 at BH q<.05 for drug vs DMSO, combination vs single drugs and interaction (AB − A − B, full-dose additive) |
| sets that survive BH | keratan sulfate biosynthesis (interaction, 8 h, q = 0.027); degradation of s2l2fn2m2masn and keratan sulfate degradation (up, PI+PD vs PI at 24 h q = 0.03; PI+5Z vs PI at 8 h q = 0.03–0.05; PI+PD vs DMSO 24 h q = 0.03); UDP-GalNAc and N-acetylglucosamine synthesis (PI+5Z interaction, 24 h, q = 0.03) |
| spermidine synthesis, PI+PD vs PI at 24 h | z = −3.4 (more down) but q = 0.41 |
| cell-cycle signature (97 genes), mean log2FC vs DMSO at 24 h | PI +0.63, PD −0.85, 5Z +0.07, PI+PD −0.84, PI+5Z +0.26 (growth at 48 h: PI 52–73%, PD 47%, 5Z 59%, PI+PD 78%, PI+5Z 92% inhibited) |
Reading: the dose-matched, replicated design does not reveal robust task-level metabolic synergy either (at most a few co-regulated glycan/keratan-sulfate sets); the one signal that appears in both datasets and in the paper's TIDE-essential "new synergy" list is up-regulated glycan / keratan sulfate degradation for PI3Ki+TAKi and PI3Ki+MEKi. The cell-cycle signature again does not follow the growth ranking.
Caveats: 3 replicates, low power; values are normalised log2 expression, not raw counts; per-condition dosing was inferred, not read from metadata; sample 'PD_1_Rep1' and 'untreated_0_Rep1' are absent from the file.

## 9. Phosphoproteomics of the same system (Tsirvouli repo; `phospho_analysis.py`)
Files: `All_ratios_Limma.txt` (4.6 MB), `DPPs_long_fdr005.csv` (2.9 MB) in `tsirvouli/`. 18,491 phosphosites on 3,927 genes, 560 Human-GEM enzymes (11.3% of sites); 0.5 / 2 / 8 h vs DMSO for PI, PD, 5Z, PI+PD, PI+5Z.
| check | result |
|---|---|
| significant sites (adj.P<.05, |logFC|≥.5) that sit on enzymes | 9.8–11.2% vs 11.3% of all sites; odds ratio 0.84–1.02 → no enrichment on metabolic enzymes |
| phospho vs mRNA on enzymes (Spearman, same condition) | −0.13 … +0.12 at every time pair → independent layers |
| PI+PD 8 h: enzymes with a significant phosphosite change (|logFC|≥.5) | 227; 37% of them have |mRNA change (24 h)| < 0.25; with same-time mRNA (8 h, VSN-scale difference) 44% (this is a small-change fraction, not a test of 'no change') |
| task level: mean enzyme phospho logFC (8 h) vs mRNA task score (24 h) | rho −0.25 … −0.01 across ≈ 88 tasks; cholesterol synthesis: phospho −0.78, mRNA +0.96 |
| sites significant in the combination where it is more extreme than BOTH single drugs | PI+PD 30 min 9%, 2 h 16%, 8 h 30% (enzymes 33%); PI+5Z 21%, 30%, 31% (enzymes 21–33%) |
| key enzymes | PHGDH, PSPH, SHMT2, ASNS, ODC1, SRM, AMD1, TXNRD1: no phosphosite measured; PSAT1 S331 +0.9 to +1.3 in PI, PD, PI+PD |
Reading: signalling-level combination effects are stronger than single-drug effects for about a third of sites, but enzyme phosphorylation is not enriched and does not track mRNA; the TIDE scores therefore describe one layer only.
Caveat: thousands of sites pass adj.P<.05 (global shifts); sites are not normalised to protein abundance; most highlighted enzymes have no measured site.

## 10. External evidence (`depmap_stream.py`, `external_check.py`)
DepMap 24Q4 CRISPR (CRISPRGeneEffect, figshare mirror of the public release; 1,178 lines, AGS = ACH-000880) and CCLE metabolomics (Li et al. 2019; 928 lines × 225 metabolites; `data.broadinstitute.org/ccle`), saved summaries in `external/`.
| gene (Chronos in AGS; < −0.5 = dependent) | value | gene | value |
|---|---|---|---|
| ASNS | −0.06 | RRM1 / RRM2 | −2.28 / −2.00 |
| PHGDH / PSAT1 / PSPH | −0.01 / −0.25 / 0.00 | DUT | −2.20 |
| ODC1 / SRM / AMD1 | −0.15 / −0.10 / −0.41 | DHODH / UMPS / CAD | −1.54 / −1.65 / −1.07 |
| TXNRD1 | −0.21 | IMPDH2 / MAT2A | −1.17 / −1.77 |
| COQ3 / COQ2 / COQ5 | −0.50 / −0.53 / −0.51 | TYMS / PPAT / PAICS | −0.29 / −0.62 / −0.77 |
- Task-essential genes (the TIDE-essential sets) are dependent in AGS more often than other metabolic genes: 29.2% vs 12.7% (AUC about 0.63–0.65) → the essential-gene idea carries real information.
- Genes of the authors' PI3Ki+MEKi "synergy-called" tasks are much more essential (median of the minimum Chronos −1.63 vs −0.29, p<0.001); the most essential tasks are the dNTP synthesis tasks (Chronos −2.3).
- mRNA of nucleotide-synthesis genes (mean log2FC vs DMSO, independent time course): de novo purine at 8 h PI −0.41, PD −0.50, 5Z −0.51, PI+5Z −0.67, PI+PD −0.74; pyrimidine −0.51 / −0.53 / −0.48 / −0.91 / −1.00; at 24 h the PI-containing conditions rebound (+0.05 to +1.2 for dNTP / RNR-TYMS) while PD and PI+PD stay down. Exploratory, not tested for significance.
- Serine/glycine synthesis genes fall in both datasets (24 h: PI −2.3 / combinations −2.4 to −3.0 in GSE285616; PI −0.8 / −1.1 to −2.0 in the time course).
- CCLE baseline for AGS (percentile among all lines): serine 19th, methionine 23rd, asparagine 51st, glycine 71st, glutamine 82nd, aspartate 85th, ornithine 96th.
Reading: the metabolic nodes the paper emphasises (serine/glycine, polyamine, thioredoxin, asparagine) show no strong baseline dependency in AGS in DepMap; CoQ synthesis genes sit at the −0.5 cutoff (COQ2 −0.53, COQ3 −0.50, COQ5 −0.51; COQ4 −0.95, PDSS1 −0.55), i.e. borderline-to-dependent, not 'non-essential'; the growth-relevant metabolic axis is nucleotide synthesis, which is suppressed early and (for MEK-containing conditions) persistently.
Caveats: DepMap screens and CCLE profiling used other media than Ham's F12 (nutrient uptake can mask dependency); metabolomics is baseline, not drug-treated; gene-set mRNA means above were not tested.

## 11. Review round 1 (external reviewer: GPT) — what I checked and what changed
Each point was verified against my files or the cited sources where possible.
| point | check | verdict / change |
|---|---|---|
| CAMERA df: estimated ρ should use min(residual df, G−2) | limma `camera.default` source: `df.camera <- min(df.residual,G-2L)` for estimated correlation, `G-2L` for fixed | reviewer correct. Re-run with df = 15: dose-scaled PI3Ki+MEKi interaction p<.05 40 → 33 tasks, BH 0 (unchanged); PI3Ki+TAKi 11 → 9, BH 0; AB vs PI3Ki alone 8 (all positive) / 3; AB vs MEKi alone 93 → 85 (BH 87 → 46). Conclusion unchanged, numbers in sections 4–7 (df≈G−2) are slightly optimistic |
| "TIDE 41 → 0 after calibration" is not like-for-like | task sets: 187 tasks, 56 have <5 genes, 11 have >200 → 120 tested, 116 distinct gene sets; gene-union, no GPR logic | reviewer correct. Correct wording: "in a separate 120-task gene-union competitive test no task passes BH 5%" |
| null/spike-in calibration covers the final pipeline | the FPR/power numbers use unmasked OLS t→z; real-DESeq2 null was degenerate (0 DEGs) | reviewer correct; keep as "OLS synthetic null" |
| 0.2–0.6 "resolution" | computed as quantiles of random-gene-set masked LFC means, no replicate-based MDE | reviewer correct. Relabel: "random-gene null quantile; power for small effects not assessed" |
| TIDE(AB−A) ≠ TIDE(AB)−TIDE(A) | min/absmax GPR and per-contrast masking are nonlinear; my contrast algebra check was gene-level only | reviewer correct. HSA/interaction on task scores are contrast-then-project quantities, not score differences |
| w-scan is a Bliss/Loewe expectation | it multiplies mRNA contrasts by a growth-derived w | reviewer correct: sensitivity analysis only. BH 1 survives at some w (w 0.6–0.8), so "disappears" was overstated: "sensitive to w" |
| "transcriptome does not proxy growth" | compared only the S+G2/M signature with digitized growth values; no task-score prediction test; different experiments | reviewer correct. Wording: "linkage not verified"; growth synergy itself rests on Flobak 2015 (CI analysis) |
| dose in Tsirvouli | Methods give the concentrations (above) | reviewer correct; slope diagnostic only used for GSE285616 (about half-dose) |
| Tsirvouli matrix = VSN | official script: `vsnMatrix` → `predict`, written as `gene_counts_flt_norm_wCTRL.csv` | reviewer correct; all "log2FC" there are VSN-scale differences (sections 8–10 relabelled). Official DE output `Transcriptomics/Output/DEGs_perContrast_all.csv` (44 MB) downloaded afterwards, see section 13 |
| CoQ non-essential | COQ2 −0.530, COQ3 −0.502, COQ5 −0.514 pass my own −0.5 cutoff; COQ4 −0.949, PDSS1 −0.548 | reviewer correct; section 10 corrected |
| 23/48 = synergy tests passed | it is combination-vs-control BH; strong 8 q values are the same test | reviewer correct; sections 2 and "Combined reading" relabelled |
| 37% is a time mismatch | same-time (8 h mRNA) fraction = 44% | reviewer correct; both reported |
| authors' BH counts "upper bound" | BH counts depend on statistic and background; not guaranteed monotone | softened to "likely optimistic" |
| reproduction gap 59 vs 2 | authors' PI3Ki: p<.025 78, BH 59; scores r = 0.95 with mine; implied null SD (|score| / z of tail p) authors 0.19 vs mine 0.26; 10,000 permutations on both sides (authors' p have 4 decimals, 2 exact zeros) | explained in section 21 (size / zero fraction of the permuted vector); the DE step and the scores were not the cause |
Null universe matters (10,000 permutations, p<.025 → BH q<.05, my reproduction):
| condition | all-genes null | model-gene null (MTEApy default) | authors' own p |
|---|---|---|---|
| PI3Ki | 44 → 2 | 37 → 2 | 78 → 59 |
| MEKi | 16 → 0 | 29 → 3 | 19 → 0 |
| TAKi | 22 → 0 | 14 → 0 | 31 → 0 |
| PI3Ki+TAKi | 29 → 3 | 40 → 5 | 56 → 3 |
| PI3Ki+MEKi | 34 → 3 | 48 → 20 | 64 → 27 |
So the PI3Ki+MEKi gap (20 vs 27) is much smaller than the PI3Ki gap (2 vs 59); section 1's headline depends on the null choice.
Done after permission (sections 12–13): official limma logFC for the Tsirvouli checks; limma-voom CAMERA / FRY in R 4.6.1. Still not possible: same-time growth, half-dose single-drug RNA (not public).

## 12. limma-voom CAMERA and FRY in R (`r_camera.R`, `r_wscan.R`, `r_null.R`; R 4.6.1, limma 3.68.5, edgeR 4.10.5)
Same counts (10,227 genes with ≥10 reads), design `~ batch + condition`, TMM + voom, the same 120 task gene sets (116 distinct), same contrasts. `camera` with `inter.gene.cor` = 0.01 (limma default, "fixed") or NA (estimated per set); `fry` = self-contained rotation test (H0: no gene of the set is differential). BH over the 120 sets.
| contrast | camera, ρ=0.01: p<.05 / BH | camera, ρ estimated: p<.05 / BH | fry: p<.05 / BH | my CAMERA-like (df=15), p<.05 / BH |
|---|---|---|---|---|
| PI3Ki vs control | 55 / 31 | 4 / 0 | 107 / 105 | 11 / 0 (df≈G−2 version) |
| PI3Ki+MEKi vs PI3Ki | 20 / 9 | 8 (all up) / 0 | 51 / 39 | 8 (all up) / 0 |
| PI3Ki+MEKi vs MEKi | 99 / 98 | 79 (all down) / 37 | 107 / 107 | 85 / 46 |
| PI3Ki+MEKi dose-scaled interaction (w=.5) | 83 (80 down) / 81 | 30 (28 down) / 0 | 90 / 84 | 33 (31 down) / 0 |
| PI3Ki+TAKi vs PI3Ki | 42 / 28 (1 down) | 4 / 0 | 33 / 18 | 3 / 0 |
| PI3Ki+TAKi dose-scaled interaction (w=.5) | 71 (68 down) / 63 | 5 / 0 | 31 / 4 | 9 / 0 |
Calibration of the real tools on this dataset (40 label shuffles within batch, 120 sets; the shuffle leaves real treatment effects in the residuals, which raises residual correlation and so overstates the fixed-ρ error somewhat):
| method | false-positive rate at p<.05 | mean BH rejections per null contrast |
|---|---|---|
| camera, ρ = 0.01 (default) | 0.28–0.35 | 24–29 |
| camera, ρ estimated | 0.02–0.03 | 0 |
| fry | 0.05–0.07 | 1.4–4.2 |
Weight scan, interaction AB − w·A − w·B (BH<.05 tasks; camera ρ=0.01 / camera ρ estimated / fry):
| w | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|
| PI3Ki+MEKi | 81 / 0 / 84 | 81 / 0 / 65 | 77 / 0 / 36 | 64 / 0 / 33 | 3 / 0 / 25 | 5 / 0 / 19 |
| PI3Ki+TAKi | 63 / 0 / 4 | 41 / 0 / 1 | 33 / 0 / 2 | 13 / 0 / 3 | 11 / 0 / 6 | 2 / 0 / 14 |
Authors' PI3Ki+MEKi synergy calls (40 testable): 39 are BH-significant under camera ρ=0.01 (but 81 of 120 sets are), 0 under camera ρ estimated, 40 under fry.
Reading:
- The estimated-ρ limma camera reproduces my CAMERA-like counts and directions closely (e.g. dose-scaled PI3Ki+MEKi 30 vs 33, AB vs PI3Ki 8 vs 8 all up, AB vs MEKi 79 vs 85), which supports the earlier calibration.
- limma's default ρ = 0.01 is anti-conservative on this dataset (FPR about 0.3; real-data estimated median ρ ≈ 0.10); the 81-task default-camera "interaction" is not reliable. Whoever applies default camera here will see strong "synergy".
- fry (nearly calibrated) shows that for PI3Ki+MEKi many task gene sets deviate from the half-dose additive expectation (33–36 sets BH at w = 0.7–0.8, 84 at w = 0.5), but the competitive test with estimated correlation finds no set that deviates MORE than other genes (BH 0 at every w). The deviation is therefore broad (a global transcriptional non-additivity), not task-specific. For PI3Ki+TAKi fry shows little at w = 0.5–0.8 (1–4 sets).
- FRY and CAMERA answer different questions; neither tests synergy in growth.

## 13. Tsirvouli checks repeated with the official limma-voom logFC (`tsirvouli_official_check.py`, `DEGs_perContrast_all.csv`)
Scale comparison (13,554 shared genes): official logFC vs VSN-scale difference r = 0.82–0.90; VSN-scale values are compressed: slope 0.50–0.64 over all genes, 0.26–0.37 for genes with mean < 7, 0.86–0.90 for mean ≥ 7.
| check | official logFC | VSN difference (earlier) |
|---|---|---|
| 24 h replication vs GSE285616 (r) | PI 0.61, PD 0.86, 5Z 0.54, PI+PD 0.84, PI+5Z 0.77 | 0.60, 0.79, 0.54, 0.77, 0.72 |
| dose slope PI+PD ~ PI / ~ PD, 8 h | 0.73 / 0.72 | 0.94 / 0.76 |
| dose slope PI+5Z ~ PI / ~ 5Z, 8 h | 0.58 / 0.83 | 0.80 / 0.95 |
| 24 h slopes (PI+PD~PI, PI+PD~PD, PI+5Z~PI, PI+5Z~5Z) | 0.46, 0.86, 0.80, 0.59 | 0.49, 0.88, 0.88, 0.77 |
| ASNS 24 h (PI, PD, 5Z, PI+PD, PI+5Z) | −1.90, −1.18, −0.22, −3.83, −2.79 | −1.83, −1.15, −0.23, −3.78, −2.74 |
| cell-cycle signature 24 h (PI, PD, 5Z, PI+PD, PI+5Z) | +0.63, −0.94, +0.05, −0.94, +0.26 | +0.63, −0.85, +0.07, −0.84, +0.26 |
| nucleotide groups 8 h, PI+PD (purine / pyrimidine / dNTP) | −0.83 / −1.09 / −0.58 | −0.76 / −1.00 / −0.52 |
| serine–glycine group 24 h (PI, PI+PD, PI+5Z) | −1.23, −2.75, −1.78 | −1.18, −2.69, −1.74 |
Reading: the qualitative statements of sections 8–10 are unchanged; VSN-scale differences understate magnitudes by 0–15% for well-expressed genes and more for low-expressed genes. With the official logFC the dose slopes in the full-dose dataset are 0.46–0.86 at 2–24 h (0.16–0.58 at 1 h), the same range as GSE285616 at 24 h (0.47–0.75). The earlier statement that the GSE285616 slopes indicate half dose is therefore NOT supported by this diagnostic; it cannot discriminate dose (selection on a noisy x biases slopes down). The half-dose design rests on the paper's Methods and Flobak 2015, not on this diagnostic. Gene-level set tests of section 8 were run on VSN data with limma-type moderated t (valid on that scale); they were not re-run on voom.

## 14. Is MTEApy/TIDE trustworthy? Tests and application to pancreatic cancer (PDAC)
Data (all public): GSE238198 (4 PDAC lines CAPAN2, CFPAC1, MiaPaCa-2, PANC-1; trametinib 3/12/24 h; 3–4 replicates; size-factor-normalised counts; limma-trend on log2(x+1), genes with mean normalised count > 5); GSE149566 (MiaPaCa-2, PANC-1 treated with fluvastatin / lovastatin / simvastatin at 2 and 20 µM, 3 replicates; log2 FPKM); DepMap 24Q4 CRISPR, 47 pancreatic lines (`external/crispr_pancreas_rows.tsv`). TIDE as in section 1 (mask padj ≥ 0.05, model-gene null, absmax, 2,000 / 1,000 permutations). Code: `r_pdac_de.R`, `pdac_tide.py`, `tide_trust_a.py`, `tide_trust_b.py`, `r_statin_de.R`, `statin_tide.py`, `statin_eval.py`.
### 14.1 Tool checks
| check | result | verdict |
|---|---|---|
| positive control: statins (HMGCR inhibitors) should raise cholesterol / mevalonate / farnesyl-PP synthesis tasks (X138, X140, X141; fixed before looking; SREBP2 feedback) | 20 µM: PANC-1, all three statins, ranks 1–4 of 187 (p ≤ 0.009, most 0.000); MiaPaCa-2: mevalonate synthesis rank 2–3, cholesterol 30–45, farnesyl-PP 11–47. 2 µM: PANC-1 cholesterol rank 6–11 (p 0.035–0.064), mevalonate score 0 (no significant genes); MiaPaCa-2 not recovered | works for a strong canonical response, weak/absent for subtle responses |
| direction | score goes UP although the enzyme is inhibited | TIDE reports the mRNA response (compensation), not the direction of flux |
| replication, independent AGS experiments (task-level Spearman, 187 tasks) | matched drug 0.60–0.72 (mean 0.66); unmatched drugs mean 0.53 (max 0.69); TIDE-essential matched 0.56 (min 0.34) vs unmatched 0.44 (max 0.76) | reproducible, but weakly drug-specific; TIDE-essential much weaker |
| added value of GPR/reaction projection | TIDE vs naive mean of masked logFC of the task genes: Pearson 0.88–0.93, top-15 overlap 7–12, no sign conflicts among tasks with both |score| > 0.2 | the network structure adds little beyond averaging the genes |
| p-values (gene-label permutation) | label-shuffle false-positive rate about 0.39 for task sets at p<.05 (section 6) | not usable as significance; use score and rank |
| low-signal inputs | PDAC 3 h: 0–50 non-zero model genes; PANC-1 24 h: 192 | with few significant genes the score is close to uninformative |
| replication across cell lines (PDAC trametinib, 24 h) | mean pairwise Spearman 0.20 (CAPAN2–CFPAC1 0.61; PANC-1 ≈ 0 with the others); 12 h vs 24 h within line 0.18–0.57 | strongly cell-line dependent |
| cross-cancer drug specificity | PDAC trametinib 24 h vs AGS MEKi 0.42–0.51 (CAPAN2, CFPAC1, MiaPaCa-2; PANC-1 ≈ 0) vs AGS PI3Ki 0.20–0.34, TAKi 0.10–0.30 | MEKi-like pattern is shared across cancers |
Overall: usable as a hypothesis generator for strong transcriptional responses (in a nutrient-feedback sense), not as a measurement of flux, synergy, or significance.
### 14.2 PDAC trametinib application (hypotheses, not results)
Number of DEGs (padj < .05): 24 h CAPAN2 2,268, CFPAC1 2,532, MiaPaCa-2 2,001, PANC-1 884.
Consensus tasks at 24 h (|score| ≥ 0.15, same sign in ≥ 3 of 4 lines): 24 tasks (12 down, 12 up); 21 of 24 have the same direction in both AGS MEKi datasets.
| task | PDAC mean score | AGS MEKi (GSE285616 / Tsirvouli) | pancreatic lines dependent (task-essential genes, min Chronos ≤ −0.5) |
|---|---|---|---|
| Oxidative phosphorylation, complex I | −0.77 | −0.41 / −0.46 | 96% (median −0.89) |
| dCTP / dTTP / dATP synthesis | −0.33 / −0.24 / −0.22 | −0.33 / −0.05 / −0.67 and −0.44 / −0.21 / −0.54 | 100% (median −3.05) |
| CTP / UTP / GTP synthesis | −0.25 / −0.21 / −0.20 | −0.29 / −0.33 / −0.62 and −0.38 / −0.38 / −0.55 | 96% / 96% / 100% |
| spermidine synthesis | −0.51 | −0.49 / −0.63 | 60% |
| cholesterol synthesis (up) | +0.64 | +0.67 / +0.85 | 11% |
| farnesyl-PP synthesis (up) | +0.59 | +0.71 / +0.81 | 89% |
Reading: MEK inhibition is associated with lower nucleotide-synthesis and complex-I transcripts and a compensatory rise of the isoprenoid arm, consistently across 3 PDAC lines and 2 AGS datasets; the nucleotide tasks and farnesyl-PP synthesis are growth-essential in pancreatic lines in DepMap.
Caveats: the S-phase proliferation signature contains TYMS, RRM1, RRM2, so part of the nucleotide signal is the proliferation arrest itself, not an independent metabolic response; PANC-1 does not follow; consensus tasks overall are not more essential than other tasks (median min-Chronos −0.33 vs −0.38); no combination, flux or metabolite data; normalised-count input and limma-trend instead of raw-count voom.

## 15. Flux-direction test: does expression-constrained FBA recover measured growth inhibition? (`eflux_growth.py`, `eflux_grid.py`)
Method: E-Flux (reaction capacity proportional to GPR expression; AND = min, OR = sum) on Human-GEM v1.18 with a Ham's-F12-like medium (glucose, 20 amino acids, vitamins, hypoxanthine, thymidine, putrescine, lipids, O2; exchange of all other metabolites closed), objective = generic human biomass minus its lipid and cofactor/vitamin pools (those pools could not be produced from the medium: vitamin-A derivatives, cobamide, phospholipid pools), GLPK, 24 h samples of GSE285616 (4 batches × 6 conditions), growth of each treated sample divided by the batch-matched control. Compared with measured 48 h growth inhibition (Flobak 2015, digitised): PI3Ki 62.5% (52–73), MEKi 47%, TAKi 59%, PI3Ki+MEKi 78%, PI3Ki+TAKi 92%.
Free parameters (not identifiable from the data): capacity scale Vmax relative to the median reaction expression, and a floor for undetected genes. 3 + 4 scale values × 3 floors scanned (18 settings).
| setting | result |
|---|---|
| zero-count genes without a floor | 600–700 reactions blocked per sample; most samples give growth 0 (dropout, not biology) |
| predicted growth inhibition, all 18 settings | between −19% and +5% (i.e. mostly predicted growth increases or no change) vs measured 47–92% |
| Spearman of predicted vs measured (n = 5) | −1.0 … +0.8, driven by differences of a few percent; best settings +0.6 to +0.8 with 7–8 of 10 concordant pairs but with the wrong sign and magnitude |
| implementation sanity: knock down RRM1+RRM2 (ribonucleotide reductase), TYMS, DHODH, PHGDH/PSAT1/PSPH, ASNS, ODC1, HMGCR, GAPDH by 90%, or 200 random metabolic genes by 50% in a control sample | growth change 0.0% in every case, also at the tightest capacity scale | growth is not limited by these enzymes in this setting; the model is insensitive to them |
Reading: this expression-constrained FBA does not recover the growth inhibition, is insensitive to the nucleotide-synthesis enzymes the earlier sections single out, and its output depends on two free scale parameters. It gives no support for inferring flux direction (or growth) from 24 h mRNA.
Limits of the test: one generic model, uncalibrated medium and uptake rates, biomass without lipid and cofactor pools, E-Flux only (not iMAT, tINIT, MOOMIN, or a context-specific reconstruction), n = 5 conditions, digitised growth values. It falsifies this naive pipeline; it does not show that no flux method can work.

## 16. Growth-level synergy with standard models (`flobak_synergy.py`; Flobak 2015 S1 Text Figs 2–38)
Data: pixel-level digitisation of the 21 representative combination growth curves (AGS, xCELLigence; value 48 h after drug addition, normalised cell index increment over the control) and of the fitted 48 h dose-response profiles of PI103, PD0325901, (5Z)-7-oxozeaenol, PKF118-310 and AKTi (CT99021 and BIRB0796 have no response in AGS). Design: single drugs at about GI50, combinations at 0.5 × of each; the authors call 4 of 21 pairs synergistic (PI3K+MEK, MEK+AKT, AKT+TAK1, PI3K+TAK1). No dose matrix exists, so MuSyC/BRAID are not applicable.
Metrics (f = fractional inhibition of the cell-index increment): HSA margin = f(AB) − max(f_A, f_B at full dose) (authors' criterion); Bliss excess = f(AB) − [1 − (1−f_A,half)(1−f_B,half)] with the half-dose single effect taken from the experiment's full-dose effect × the half/full ratio of the fitted dose-response curve; Loewe combination index CI = Σ d_i / Dx_i with Dx_i the single-drug dose giving the combination effect on its curve (a "direct" variant reads half-dose effects straight from the separate dose-response experiment).
| combination (Fig) | f_A / f_B (full dose) | f(AB) | HSA margin | Bliss expected | Bliss excess | Loewe CI |
|---|---|---|---|---|---|---|
| PI3Ki + MEKi (18) | 0.52 / 0.48 | 0.77 | +0.25 | 0.59 | **+0.18** (direct +0.14) | **0.26** (direct 0.37) |
| PI3Ki + TAKi (21) | 0.76 / 0.61 | 0.96 | +0.20 | 0.71 | **+0.25** (direct +0.38) | **0.25** (direct 0.18) |
| MEKi + AKTi (19) | 0.42 / 0.29 | 0.83 | +0.41 | 0.38 | +0.45 | 0.26 |
| AKTi + TAKi (20) | 0.25 / 0.61 | 0.89 | +0.28 | 0.46 | +0.43 | 0.33 |
| PI3Ki + AKTi (22, non) | 0.76 / 0.25 | 0.61 | −0.15 | 0.53 | +0.08 | 1.09 |
| MEKi + TAKi (26, non) | 0.47 / 0.67 | 0.55 | −0.12 | 0.67 | −0.12 | 1.05 |
Class separation over all 21 pairs (4 called synergistic, 17 not): HSA margin synergistic +0.20…+0.41 vs non −1.37…−0.10 (AUC 1.00); Bliss excess +0.18…+0.45 vs −0.42…+0.11 (AUC 1.00; direct variant 0.97); CI 0.25–0.33 vs 0.60–2.07 for the 10 non-synergistic pairs with a defined CI (all lower; CI is undefined when f(AB) < 0.1).
Reproducibility context: the same single drug at the same dose gives f = 0.52–0.76 (PI103, sd 0.10), 0.42–0.64 (PD0325901, sd 0.08), 0.41–0.67 (5Z, sd 0.11) across experiments; the Bliss excess of the 17 non-synergistic pairs is −0.13 ± 0.18 (mean ± sd), so the PI3Ki+MEKi excess is about 1.8 sd and PI3Ki+TAKi about 2.2 sd above that reference.
Sensitivity of the Bliss excess to the half/full dose effect ratio r (curves give 0.66 PI103, 0.80 PD0325901, 0.69 5Z): r = 0.5 / 0.66 / 0.8 / 0.9 / 1.0 → PI3Ki+MEKi +0.34 / +0.22 / +0.13 / +0.08 / +0.02; PI3Ki+TAKi +0.39 / +0.26 / +0.16 / +0.11 / +0.06. The excess vanishes only if halving the dose barely reduces the effect (r → 1), which the fitted curves contradict.
Reading: at the growth level, both combinations are synergistic by HSA, Bliss and Loewe criteria with moderate size (Bliss excess about 0.2; CI about 0.25), and the classification agrees with the authors' calls across the 21 pairs. PI3Ki+TAKi is stronger and more certain than PI3Ki+MEKi. This is the real synergy; it is not shown by the transcriptome-based task scores (sections 5–8).
Caveats: one representative replicate per figure; colour and axis readings from raster images (about ±0.05 cell index); single-drug half-dose effects are not measured in the combination experiments (derived from separate dose-response fits); no confidence intervals for the CI; MuSyC/BRAID need a dose matrix that does not exist here.

## 17. Development step: response-to-dependency (R2D) map for MEK inhibition in pancreatic cancer (`pdac_r2d.py`, `figs/pdac_r2d.png`)
Purpose: turn the drug-response mapping idea of the paper into a partner-target prioritisation that does not rely on TIDE scores or gene-label permutation p-values. Gene-level inputs only: trametinib 24 h (GSE238198; CAPAN2, CFPAC1, MiaPaCa-2, PANC-1; limma-trend logFC), dependency in 47 DepMap pancreatic lines (Chronos), AGS MEKi (GSE285616 and Tsirvouli, official logFC), simvastatin in PANC-1 and MiaPaCa-2 (GSE149566).
Rule: response up (mean log2FC ≥ 0.5, padj < .05 with the same sign in ≥ 3 of 4 lines) or down (≤ −0.5, ≥ 3 of 4 lines); essential = median Chronos ≤ −0.5 and ≥ 50% of pancreatic lines dependent. 1,710 metabolic genes with data.
| class | n | genes |
|---|---|---|
| A: up and essential (escape route) | 4 | HMGCS1 (median Chronos −1.95, 100% of lines), SCD (−0.99), MVK (−0.88), MVD (−0.51); HMGCR is up in 1 line only (+0.51; −1.10, 98%) |
| B: down and essential (amplified vulnerability) | 17 | FOSL1 (−2.75), RRM2 (−1.24; −2.37), RRM1 (−0.83; −3.05), CTPS1, GMPS, POLA1/2, RFC5, PPAT and other replication genes |
| essential but no response | – | FNTA (−1.05), FNTB (−0.70), GGPS1 (−1.33; prenyltransferases) |
| up but not essential | 40 | cholesterol branch: SQLE, LSS, DHCR7, FDFT1, MSMO1, FASN … |
Cross-checks: every class A gene goes up in both AGS MEKi datasets (100%) and in simvastatin-treated PDAC lines (statin logFC ≥ 0.5 for 4 of 4, vs 14% of all metabolic genes), i.e. the A genes are the SREBP2 feedback program; class B genes go down in AGS MEKi (82% / 100%) and mostly also under simvastatin (6% up).
Pattern: the essential part of the induced program is the shared mevalonate trunk (HMGCS1, HMGCR, MVK, MVD → FPP/GGPP; prenyltransferases FNTA/FNTB/GGPS1 are essential in 83–98% of lines), while the cholesterol branch is induced but dispensable in culture.
Literature context (web search, 2026-10-12, not verified in depth): statin + MEK inhibitor (CH5126766, trametinib) induced apoptosis with mevalonate/GGPP rescue in a 2018 Oncotarget study (cancer type outside PDAC); SREBP2-dependent HMGCR/HMGCS1 induction after statins is established. No source found that shows MEK-inhibitor-induced mevalonate-trunk induction in PDAC; trametinib + gemcitabine or hydroxychloroquine failed to show benefit clinically.
Hypotheses (not results): H1 trametinib induces the mevalonate trunk through SREBP2 in PDAC lines; H2 an HMGCR inhibitor (or a prenylation inhibitor) combined with trametinib is synergistic in PDAC lines and the effect is reversed by mevalonate or geranylgeraniol but not by cholesterol; H3 (amplified arm) trametinib lowers RRM1/RRM2/CTPS1 and dNTP/CTP pools; its combination with nucleotide-targeting drugs depends on sequence and dose.
Limits: logFC from log2(normalised count + 1) with limma-trend, genes with mean count < 5 removed; mRNA only; essentiality measured in standard culture (serum supplies cholesterol); thresholds chosen once, not tuned; class A has 4 genes only (the pattern is one pathway, not four independent findings); no combination or flux data.

### 17b. Which parts of Benedicto 2025 the R2D step uses, and what the paper's own data add
| part of the paper | used for |
|---|---|
| GSE285616 (AGS, 5 conditions at 24 h, incl. PI3Ki+MEKi) | the same escape-route question asked under the paper's own best combination; AGS as the discovery system, pancreatic lines as target |
| Human-GEM metabolic tasks (187) and the task-essential gene sets of TIDE-essential | annotation of A/B genes by task; essential-gene prior |
| Supplementary Data 3 (task scores) | the isoprenoid tasks are already in the paper's tables |
| the PI3K+MEK combination | the rationale: KRAS acts through RAF–MEK–ERK and PI3K |
Paper's Supplementary Data 3 (TIDE score, p; sheets as named in the file, the TAKi/MEKi header cells are swapped there): cholesterol synthesis (X138) MEKi-sheet 0.70 (p 0.006), TAKi-sheet 0.68 (p 0.004), PI3Ki 0.23 (p 0.025), PI3Ki+MEKi 0.37 (p 0.038), PI3Ki+TAKi 0.32 (p 0.020); mevalonate synthesis (X140) 0.84 / 0.94 / 0.44 / 0.60 / 0.61 (p 0.044 / 0.030 / 0.062 / 0.080 / 0.057); farnesyl-PP synthesis (X141) 0.67 / 0.70 / 0.29 / 0.21 / 0.33 (p 0.041 / 0.024 / 0.078 / 0.258 / 0.098). The up-regulation of the isoprenoid arm is present in the paper's own results.
Gene level in AGS (log2FC, * padj < .05; PI3Ki | MEKi | TAKi | PI3Ki+TAKi | PI3Ki+MEKi): HMGCS1 +0.53* | +1.48* | +1.19* | +0.64* | +0.79*; HMGCR +0.41* | +0.84* | +0.93* | +0.56* | +0.57*; MVD +0.36* | +1.60* | +1.42* | +0.51* | +0.65*; MVK +0.43* | +0.66* | +0.58* | +0.04 | +0.02; FDPS −0.02 | +0.57* | +0.68* | +0.28* | −0.01; SCD −0.35* | +0.33* | +0.20 | −0.33* | −0.47*; SREBF1 −0.80* | −0.17 | −0.21* | −0.58* | −0.71*.
Reading: under the paper's half-dose PI3Ki+MEKi combination the mevalonate trunk (HMGCS1, HMGCR, MVD, IDI1) stays induced at about half the MEKi-alone magnitude, which is what half dosing alone predicts (only MVK and FDPS are clearly lower); the SCD arm is reversed by PI3K inhibition (SCD −0.47, SREBF1 −0.71). So the escape that remains under the combination is the mevalonate trunk, not the SCD arm.
Dependency of the same genes (Chronos, AGS | gastric median | pancreatic median): HMGCS1 −1.36 | −1.81 | −1.95; HMGCR −0.87 | −0.99 | −1.10; MVK −1.04 | −0.81 | −0.88; GGPS1 −1.58 | −1.36 | −1.33; FNTA −1.25 | −1.02 | −1.05; SQLE −0.24 | −0.09 | −0.14; LSS −0.03 | −0.03 | −0.05. The essential trunk and the dispensable cholesterol branch are the same in gastric and pancreatic lines, so the AGS system carries over.
Paper's task-essential gene sets: HMGCR ∈ X140 mevalonate synthesis; MVK, MVD ∈ X141 farnesyl-PP synthesis; RRM1, RRM2 ∈ X13–X16 dNTP synthesis; GMPS ∈ X11, X15; PPAT ∈ X9, X11, X13, X15 (HMGCS1 and CTPS1 are not in the paper's essential sets). In pancreatic lines the task-essential genes are dependent more often than other metabolic genes, but weakly (13.6% vs 10.9%; AUC 0.64); X141 (3 genes, minimum Chronos −0.88, 67% dependent) ranks 20th of 95 tasks.

## 18. Limitation-resolution round for the R2D hypothesis (`pdac_validate.py`)
| limitation | test | result |
|---|---|---|
| half-dose confound in the paper's combination | Tsirvouli AGS data: single drugs and combinations at the same full doses, official logFC, 1–24 h | no induction at 1–4 h, onset at 8 h (MEKi: HMGCS1 +0.83, HMGCR +0.69, MVK +0.85, MVD +0.90), strong at 24 h. At 24 h PI3Ki alone also induces the trunk (HMGCS1 +1.75, MVD +1.16, IDI1 +1.62), MEKi +1.90 / +1.72 / +1.25, PI3Ki+MEKi +2.47 / +2.07 / +1.89, i.e. PI3K inhibition does not attenuate it and the full-dose combination is higher. INSIG1 +2.0 to +2.6, SREBF2 +0.2 to +0.7, SCAP flat. SCD arm: PI +0.45, PD +0.65, PI+PD +0.06 (neutralised in the combination). LDLR is down under MEKi (−2.09) and up under PI3Ki at 1–8 h, so the program is not a pure SREBP2 response |
| single-dataset evidence | GSE206267 (KP4, MiaPaCa-2; 100 nM trametinib, 48 h, DMEM + 10% FBS; raw counts) | partial: KP4 HMGCS1 +0.46 (p<.05), HMGCR +0.28, MVD +0.62, 17 of 24 sterol-program genes up (sign test p = 0.064; program mean +0.18, not different from random 24-gene sets, percentile 0.26); MiaPaCa-2 HMGCS1 +0.82, HMGCR −0.17, program mean −0.02 |
| generality to RAS inhibition | GSE315130 (PSN1, SW1990; RAS(ON) inhibitor RMC-6236 vs DMSO; FPKM; time and dose not stated in the record) | not seen: trunk slightly down (program mean −0.19 PSN1, −0.13 SW1990; sign test p = 0.007 and 0.023 for down); with tazemetostat PSN1 +0.21 (n.s.), SW1990 −0.68 |
| program-level evidence instead of 4 genes | 24 pre-listed SREBP2/sterol genes, same-sign counts and random-set percentile | GSE238198 24 h: 22 / 23 / 18 / 23 of 24 up in CAPAN2 / CFPAC1 / MiaPaCa-2 / PANC-1 (percentile 0.000, 0.000, 0.005, 0.000); simvastatin PANC-1 24 of 24, MiaPaCa-2 17 of 24 (0.0085) |
Reading: the induction of the mevalonate trunk by MEK inhibition is strong in one 4-line PDAC dataset at 24 h and in two AGS datasets, is only partly reproduced in an independent trametinib dataset at 48 h, and is absent (slightly reversed) with the RAS(ON) inhibitor in two lines of one dataset. The hypothesis therefore holds for early MEK-inhibitor response in some lines and cannot be extended to daraxonrasib-class drugs.
Not resolved: essentiality in standard culture only; mRNA only; no PDAC data for PI3K inhibitors or combinations; the SREBP2 attribution (LDLR argues against a pure SREBP2 response); time and dose of GSE315130; effect size in the second trametinib dataset is small.

## 19. Second limitation round: kinetics, proliferation confound, signature breadth, KRAS inhibition at 24 h (`pdac_limits2.py`, inline test on GSE279749)
| question | test | result |
|---|---|---|
| time course in the same four PDAC lines (GSE238198) | mean log2FC of the 24 sterol-program genes at 3 / 12 / 24 h | CAPAN2 +0.09 / +0.63 / +1.06; CFPAC1 −0.18 / +0.55 / +1.31; MiaPaCa-2 −0.10 / +0.04 / +0.33; PANC-1 +0.05 / +0.05 / +1.12. Onset by 12 h in the two low-grade lines, at 24 h in the others; the independent 48 h dataset (section 18) is weaker, so a rise-then-decline in some lines is possible but is not shown (different labs, doses) |
| is it only a proliferation-arrest correlate? | cell-cycle signature (97 S + G2/M genes) vs program over 12 line × time points; trunk genes are not in the cell-cycle list | Spearman −0.78 (program) and −0.86 (trunk), but time drives both. At the fixed 24 h time point the strongest cell-cycle suppression is in MiaPaCa-2 (−1.46) with the weakest induction (+0.33) and the weakest suppression in PANC-1 (−0.45) has +1.12, so induction does not scale with arrest across lines (n = 4) |
| broader, data-derived signature | 100 genes most induced by simvastatin in both PDAC lines (outside the 24 listed genes; mostly stress genes: TSC22D3, S100P, ATF3, HMOX1, LCN2, NUPR1, RHOB) | up under trametinib in CAPAN2 (+0.73), CFPAC1 (+0.65), MiaPaCa-2 (+0.30) at 24 h and in both 48 h lines (KP4 +0.51, MiaPaCa-2 +1.26), not in PANC-1 (−0.25) and down in AGS (PD −0.37, PI −0.67, PI+PD −0.56). The statin response in PDAC lines contains a stress component that MEK inhibition shares; the 48 h MiaPaCa-2 response is of this type, not of the sterol program |
| KRAS inhibition, documented time and dose (GSE279749; 24 h; MRTX1133 50 nM in AsPC-1, PANC-1, SW1990; BI-2865 100 nM in CFPAC-1; 3 replicates; FPKM, 3,909 genes) | trunk and program (18 of 24 genes present) | PANC-1: HMGCS1 +0.60, HMGCR +0.27, MVD +0.36, IDI1 +0.47, program +0.33 (14 of 18 up, random-set fraction 0.000); AsPC-1 −0.10 (0.64), SW1990 −0.17 (0.47), CFPAC-1 +0.05 (0.08). The KRAS G12D inhibitor reproduces the induction in 1 of 4 lines (PANC-1 only) at 24 h; the same line also responds to trametinib at 24 h (+1.05 trunk) |
Reading: induction of the mevalonate trunk occurs after MAPK-pathway inhibition in some PDAC lines (CAPAN2, CFPAC1, PANC-1 with trametinib; PANC-1 with MRTX1133), not in others (AsPC-1, SW1990, MiaPaCa-2 weak), and was not seen with RMC-6236 in two lines. It is line-dependent rather than a general response of KRAS-pathway inhibition.
Caveats: GSE279749 has only 3,909 genes (a filtered table) and FPKM; differences in dose/time between studies; 4 lines each.

## 20. What distinguishes inducing from non-inducing PDAC lines? (`depmap_pdac_omics.py`, `pdac_responders.py`; exploratory, n = 7)
Inducibility = mean program log2FC (24 sterol/SREBP genes) after MAPK-pathway inhibition over the datasets of sections 18–19: CAPAN2 +1.06, PANC-1 +0.73, CFPAC1 +0.68, KP4 +0.18, MiaPaCa-2 +0.16, AsPC-1 −0.10, SW1990 −0.17. DepMap 24Q4 pancreatic models (55 with expression, 47 with CRISPR; CAPAN2 has no CRISPR row).
| variable | result (Spearman with inducibility, n = 7 unless noted) |
|---|---|
| baseline mean expression of the six trunk genes (HMGCS1, HMGCR, MVK, MVD, IDI1, FDPS) | **rho = +0.96**; leave-one-line-out 0.94–1.00; only 0.3% of 3,000 random 6-gene sets reach \|rho\| ≥ 0.96 (95th percentile of \|rho\| 0.82) |
| individual baseline genes | INSIG1 +0.89, MBTPS1 +0.89, SQLE +0.86, MVK +0.86, HMGCR +0.79, MVD +0.75; SPRY2 (ERK output) −0.86 (44 genes tested in 3 panels, so 2 chance hits at p<.05 are expected) |
| dependency on the SREBP processing machinery (CRISPR, n = 6) | SCAP −0.89 (p .019), MBTPS1 −0.94 (p .005), SREBF2 −0.71, MBTPS2 −0.77: the inducing lines are the ones that depend more on SCAP / site-1 protease; HMGCS1 +0.83 (inducing lines depend less on it) |
| classical-minus-basal score, KRAS allele, TP53 / SMAD4 / CDKN2A | no relation (rho 0.14; alleles G12V / G12D / G12C spread over both groups; SMAD4 and CDKN2A loss only in AsPC-1) |
| baseline-vs-induction inside single studies (n = 4 each) | GSE238198 rho 0.4; GSE279749 rho 0.4 (uninformative) |
Reading: induction is found in lines with an already high, SCAP-dependent SREBP2-type program, not in lines with a low program; this argues against a "ceiling" explanation and for line-intrinsic pathway capacity. The same two lines that were measured with a statin and with trametinib (PANC-1 strong, MiaPaCa-2 weak) rank the same way.
Caveats: 7 lines from different studies, doses and times; inducibility is a composite; trunk genes were chosen from the response data; many variables tested; no patient data. Prediction to be tested in the lab: statin synergy with trametinib should be larger in lines with high baseline trunk expression and SCAP dependency (CAPAN2, PANC-1, CFPAC1) than in low-expression lines (SW1990, AsPC-1, MiaPaCa-2).

## 21. Why the reproduction gap (BH 59 vs 2): the size of the permuted vector (`repro_null_universe.py`, `repro_null_scan.py`, `repro_null_all.py`, `repro_variants2.py`)
Paper Methods: "Genes with non-significant adjusted p-values were masked by setting their log-FC values to zero … statistical significance … through 10,000 permutations … p < 0.025". The size and composition of the permuted gene vector are not stated.
What was ruled out (PI3Ki, my pipeline):
- DE step: my DEG counts (|LFC| ≥ 1, padj < .05) 1,342 / 2,216 / 1,571 / 1,761 / 2,911 vs the authors' Supplementary Data 1 lists 1,274 / 2,106 / 1,616 / 1,769 / 2,991 (PI3Ki / MEKi / TAKi / PI3Ki+TAKi / PI3Ki+MEKi); the paper says about 2,000 per condition.
- Scores: r = 0.95 with the authors' TIDE scores (slope 0.90).
- Masking rule: padj < .05 & |LFC| ≥ 1: p<.025 22; padj < .05 & |LFC| ≥ 0.5: 23; padj < .01: 33; unshrunk & |LFC| ≥ 1: 23; baseline 34 (BH 2 in all) — none approaches 78.
- Shrinkage or null universe choices tested in section 11 — no.
Decisive test: apply the AUTHORS' scores to null distributions built from my masked log2FC vector, zero-padded to size N (gene-label permutation, 2,000 permutations):
| N (entries in the permuted vector) | non-zero fraction | tasks p < .025 (authors 78) | BH q < .05 (authors 59) | RMSE of −log10 p vs authors |
|---|---|---|---|---|
| 2,443 (model genes with data, my default) | 0.51 | 44 | 3 | 0.53 |
| 3,000 | 0.41 | 56 | 4 | 0.42 |
| 4,000 | 0.31 | 71 | 5 | 0.25 |
| **5,000** | **0.25** | **78** | **53** | **0.19** |
| 6,000 | 0.21 | 84 | 61 | 0.26 |
| 8,000 | 0.16 | 88 | 69 | 0.39 |
| 16,932 (all tested genes) | 0.07 | 106 | 94 | 0.75 |
| 23,373 (all genes in the count file) | 0.05 | 113 | 106 | – |
Other conditions (authors' scores; p<.025 / BH): PI3Ki+TAKi authors 56 / 3, best N = 4,000 gives 55 / 8 (RMSE 0.14); PI3Ki+MEKi authors 64 / 27, model-gene vector gives 55 / 27 and N = 4,000 gives 63 / 39. MEKi and TAKi (authors 19 / 0 and 31 / 0) are not matched by any N (RMSE 0.32–0.46), consistent with the weak score agreement for MEKi (r = 0.40).
Reading:
1. With the same scores, the number of "significant" tasks varies 44 → 117 and the BH survivors 3 → 107 only by changing how many zero entries sit in the permuted vector. The authors' counts are reproduced when the vector has about 4,000–5,000 entries (about 25–30% non-zero), which is neither my default (model genes with data, 51% non-zero) nor the whole transcriptome.
2. The gap is therefore a consequence of an unreported analysis choice in the permutation null, not of the DE step or of the scores. The paper's task-level significance, and the 59 PI3Ki "BH survivors" that follow from it, depend on that choice.
3. My earlier "all-genes" run (section 1) is not the same as zero-padding: it keeps the non-model genes' own non-zero log2FC (about 19% non-zero in total), which widens the null again.
Caveats: the authors' gene universe is inferred, not known; N = 4,000–5,000 matches three of five conditions; MEKi / TAKi remain unexplained; 2,000 permutations per null.

## 22. The MEKi (r = 0.40) and TAKi mismatches (`tide_from_authors_deg.py`, `repro_idmap.py`, `repro_idmap_null.py`, `repro_deg_null.py`, `taki_diff.py`)
Sample-to-condition assignment of GSE285616 was verified from the series matrix (column order: control, PI103, PD0325901, (5Z)-7-oxozeaenol, PI103+5Z, PI103+PD; batches 1–4), so the mismatch is not a label swap.
| finding | evidence |
|---|---|
| The Pearson r = 0.40 for MEKi is driven by 5 tasks | tryptophan degradation tasks X124–X128 (kynurenine pathway): authors −0.4 to −1.2, mine +1.3 to +3.9; without these 5 tasks r = 0.90 (Spearman over all tasks was already 0.75) |
| Cause: symbol-based gene mapping | in the authors' Supplementary Data 1 the symbol TDO2 carries ENSG00000262635 (log2FC +5.5, padj .006), which is not a Human-GEM gene, whereas the model's TDO2 is ENSG00000151790. My symbol join gave TDO2's +6.5 (raw counts 0,0,0,0 in control vs 5,1,9,43 in MEKi) to the model gene, so absmax in the IDO1/IDO2/TDO2 OR-rule picked TDO2 instead of IDO1 (−2.36 → task score −1.18 = −2.36 / 2 reactions). 26 symbols in the authors' DEG lists carry an ID different from Human-GEM's |
| Fix: join by the authors' Ensembl IDs (where given) | r(scores) PI3Ki 0.95 → 0.95; **MEKi 0.40 → 0.89**; TAKi 0.82 → 0.82; PI3Ki+TAKi 0.94 → 0.95; **PI3Ki+MEKi 0.78 → 0.92** |
| Scoring directly from the authors' own DEG lists (|LFC| ≥ 1, padj < .05, others 0) | r with their scores 0.95 / 0.92 / 0.84 / 0.93 / 0.94 (PI3Ki / MEKi / TAKi / PI3Ki+TAKi / PI3Ki+MEKi); the DEG-only masking is close to, but does not exactly equal, their scores (37 of 187 exact) |
| TAKi residual (r = 0.82; 0.935 without the 5 largest-difference tasks) | ultra-low-count genes: PAPSS2 (6 reads in total, a DEG +2.05 in the authors' list, removed by my count filter), IDO1 (padj 0.16 here vs 0.004 there), HGD / HPD / CAPN3 / TXNRD3 (10–40 reads). R DESeq2 + apeglm and my pydeseq2 differ for such genes; the authors' exact inputs are not public |
p-values with the corrected mapping (authors' scores against nulls, counts p<.025 / BH): MEKi authors 19 / 0, model-gene universe 30 / 0, N = 5,000 27 / 0; TAKi authors 31 / 0, 26 / 2, 38 / 1. With a null built from the authors' own DEG vectors (zero-padded; N = 3,041 = all Human-GEM genes): MEKi 20 / 0 (authors 19 / 0), TAKi 25 / 1 (31 / 0), PI3Ki+MEKi 53 / 27 (64 / 27), PI3Ki+TAKi 42 / 4 (56 / 3), PI3Ki 47 / 3 (78 / 59; 77 / 31 at N = 5,000).
Reading: the MEKi score mismatch was my gene-ID mapping, not the authors' analysis; it is fixed. The TAKi mismatch is confined to a handful of ultra-low-count genes and cannot be closed without their exact count inputs. For the p-values, MEKi and PI3Ki+MEKi are reproduced in count by a model-gene-sized null built from the authors' DEG vectors, while PI3Ki and PI3Ki+TAKi need a wider vector (about 5,000); no single universe reproduces all five conditions, i.e. the undocumented permutation universe and masking rule determine the reported counts.
Effect on earlier statements: the "r = 0.40 for MEKi" and "r = 0.78 for PI3Ki+MEKi" in section 2 are artefacts of symbol mapping; the corrected values are 0.89 and 0.92. Gene-level analyses elsewhere in this report join by symbol and are unaffected.

