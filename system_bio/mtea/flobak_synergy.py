"""Growth-level synergy from Flobak 2015 S1 Text: pixel-digitised 48 h growth (21 combinations) and single-drug dose-response profiles.
For each combination (each drug at 0.5 x its GI50-like dose): authors' HSA-type criterion, Bliss excess at the half dose, Loewe combination index (CI)."""
import numpy as np, pandas as pd
from PIL import Image
from scipy import stats
S = "si/"
def load(f): return np.array(Image.open(S + f).convert("RGB")).astype(int)
COL = {"red": lambda p: (p[..., 0] > 200) & (p[..., 1] < 70) & (p[..., 2] < 70),
       "green": lambda p: (p[..., 0] < 90) & (p[..., 1] > 190) & (p[..., 2] < 90),
       "blue": lambda p: (p[..., 0] < 70) & (p[..., 1] < 70) & (p[..., 2] > 190),
       "magenta": lambda p: (p[..., 0] > 190) & (p[..., 1] < 70) & (p[..., 2] > 190)}

# ---------------- dose-response profiles (48 h cell index vs log10 dose; fitted curves) ----------------
DR = {  # drug: (file, xmin, xmax, ymin, ymax, log10 dose used in combinations)
    "PI103": ("fl_p12_img1.png", -8.5, -4.0, 0.0, 4.8, np.log10(0.7e-6)),
    "PD0325901": ("fl_p12_img5.png", -10.0, -4.0, 1.2, 4.8, np.log10(35e-9)),
    "5Z": ("fl_p13_img1.png", -7.5, -5.0, 1.0, 3.5, np.log10(0.5e-6)),
    "PKF118-310": ("fl_p12_img3.png", -8.0, -6.0, 0.0, 5.4, np.log10(150e-9)),
    "AKTi": ("fl_p13_img3.png", -6.0, -4.0, -0.7, 3.5, np.log10(10e-6)),
}
def frame(im):
    dark = (im.sum(axis=2) < 640) & (abs(im[:, :, 0] - im[:, :, 1]) < 15)      # gray / black grid and frame (incl. anti-aliased)
    rows = dark.sum(1); cols = dark.sum(0)
    r = np.where(rows > 0.8 * im.shape[1])[0]; c = np.where(cols > 0.8 * im.shape[0])[0]
    if len(r) == 0 or len(c) == 0:
        r = np.where(rows > 0.7 * rows.max())[0]; c = np.where(cols > 0.7 * cols.max())[0]
    return c.min(), c.max(), r.min(), r.max()
def dose_curve(drug):
    f, xmin, xmax, ymin, ymax, dose = DR[drug]
    im = load(f); x0, x1, y0, y1 = frame(im)
    m = (im[..., 0] > 180) & (im[..., 1] < 170) & (im[..., 2] < 170) & ((im[..., 0] - im[..., 1]) > 80)
    xs, vals = [], []
    for x in range(x0 + 2, x1 - 1):
        ys = np.where(m[y0 + 2:y1 - 1, x])[0]
        if len(ys): xs.append(x); vals.append(ys.mean() + y0 + 2)
    vals = pd.Series(vals).rolling(9, center=True, min_periods=1).median().values
    logd = xmin + (np.array(xs) - x0) / (x1 - x0) * (xmax - xmin)
    ci = ymax - (vals - y0) / (y1 - y0) * (ymax - ymin)
    return logd, ci, dose
CURVES = {}
for d in DR:
    logd, ci, dose = dose_curve(d)
    ci0 = ci[:5].mean()
    fa = (ci0 - ci) / (ci0 - 1.0)
    fa = np.maximum.accumulate(np.clip(fa, -0.2, 1.5))     # enforce monotone increasing with dose
    CURVES[d] = (logd, fa, dose, ci0)
    print(f"{d:11} log10 dose range {logd.min():.2f}..{logd.max():.2f}, CI0 {ci0:.2f}; fa at dose {np.interp(dose, logd, fa):.2f}, at dose/2 {np.interp(dose-np.log10(2), logd, fa):.2f}, at dose/4 {np.interp(dose-np.log10(4), logd, fa):.2f}; max fa {fa.max():.2f}")
def fa_at(drug, logdose): return float(np.interp(logdose, CURVES[drug][0], CURVES[drug][1]))
def dose_for(drug, e):
    logd, fa, _, _ = CURVES[drug]
    if e > fa.max(): return np.inf
    if e <= fa[0]: return 0.0
    return 10 ** float(np.interp(e, fa, logd))        # fa monotone

# ---------------- time-course figures: value 48 h after drug addition ----------------
FIGS = {  # figure: (file, {drug: colour}, authors' call)
    18: ("fl_p18_img0.png", {"PD0325901": "green", "PI103": "blue"}, "syn"),
    19: ("fl_p18_img1.png", {"PD0325901": "blue", "AKTi": "green"}, "syn"),
    20: ("fl_p18_img2.png", {"AKTi": "green", "5Z": "blue"}, "syn"),
    21: ("fl_p18_img3.png", {"PI103": "green", "5Z": "blue"}, "syn"),
    22: ("fl_p18_img4.png", {"AKTi": "blue", "PI103": "green"}, "non"),
    23: ("fl_p18_img5.png", {"AKTi": "blue", "BIRB0796": "green"}, "non"),
    24: ("fl_p19_img0.png", {"AKTi": "green", "CT99021": "blue"}, "non"),
    25: ("fl_p19_img1.png", {"AKTi": "green", "PKF118-310": "blue"}, "non"),
    26: ("fl_p19_img2.png", {"PD0325901": "blue", "5Z": "green"}, "non"),
    27: ("fl_p19_img3.png", {"PD0325901": "green", "PKF118-310": "blue"}, "non"),
    28: ("fl_p19_img4.png", {"PD0325901": "green", "CT99021": "blue"}, "non"),
    29: ("fl_p19_img5.png", {"PD0325901": "blue", "BIRB0796": "green"}, "non"),
    30: ("fl_p20_img0.png", {"PI103": "green", "PKF118-310": "blue"}, "non"),
    31: ("fl_p20_img1.png", {"PI103": "blue", "CT99021": "green"}, "non"),
    32: ("fl_p20_img2.png", {"PI103": "blue", "BIRB0796": "green"}, "non"),
    33: ("fl_p20_img3.png", {"CT99021": "green", "PKF118-310": "blue"}, "non"),
    34: ("fl_p20_img4.png", {"CT99021": "blue", "BIRB0796": "green"}, "non"),
    35: ("fl_p20_img5.png", {"5Z": "blue", "PKF118-310": "green"}, "non"),
    36: ("fl_p21_img0.png", {"5Z": "blue", "CT99021": "green"}, "non"),
    37: ("fl_p21_img1.png", {"5Z": "blue", "BIRB0796": "green"}, "non"),
    38: ("fl_p21_img2.png", {"PKF118-310": "blue", "BIRB0796": "green"}, "non"),
}
PXH = 11.3; X0 = 53.0; T_END = 69.0; T_DRUG = 21.0
def read_fig(fn, cols):
    im = load(fn)
    xr = lambda t: int(round(X0 + PXH * t))
    def row_at(color, t, hw=2):
        m = COL[color](im)[:, xr(t) - hw: xr(t) + hw + 1]
        ys = np.where(m.any(axis=1))[0]
        return float(np.median(ys)) if len(ys) else np.nan
    r21 = np.nanmedian([row_at(c, T_DRUG - 0.8) for c in ("red", "green", "blue", "magenta")])
    scale = (460.0 - r21) / 1.5          # px per cell-index unit (bottom frame = -0.5, drug addition = 1.0)
    val = lambda row: 1.0 + (r21 - row) / scale
    out = {"control": val(row_at("red", T_END)), "combo": val(row_at("magenta", T_END))}
    for d, c in cols.items(): out[d] = val(row_at(c, T_END))
    return out, scale

rows = []
for fig, (fn, cols, call) in FIGS.items():
    v, scale = read_fig(fn, cols)
    inc = lambda x: x - 1.0
    ctrl = inc(v["control"])
    drugs = list(cols)
    f = {d: 1 - inc(v[d]) / ctrl for d in drugs}; fab = 1 - inc(v["combo"]) / ctrl
    rows.append(dict(fig=fig, A=drugs[0], B=drugs[1], call=call, control=ctrl, fA=f[drugs[0]], fB=f[drugs[1]], fAB=fab, scale=scale))
df = pd.DataFrame(rows).set_index("fig")
LH = np.log10(2.0)
def metrics(r):
    out = {}
    hA = hB = None
    for key, d, fx in (("A", r.A, r.fA), ("B", r.B, r.fB)):
        if d in CURVES:
            logd, fa, dose, ci0 = CURVES[d]
            full_c = fa_at(d, dose); half_c = fa_at(d, dose - LH)
            ratio = half_c / full_c if full_c > 0.05 else 0.5
            h_cal = max(0.0, fx) * ratio                  # within-experiment potency x curve shape
            h_dir = half_c                                # taken directly from the separate dose-response experiment
            # Loewe: dose-equivalent shift so that the curve reproduces the within-experiment full-dose effect
            dx = dose_for(d, r.fAB); s_dose = dose_for(d, max(fx, fa[0] + 1e-3))
            eq_half_cal = 0.5 * (s_dose if np.isfinite(s_dose) and s_dose > 0 else 10 ** dose)
            term_cal = 0.0 if (np.isinf(dx) or dx <= 0) else eq_half_cal / dx
            term_dir = 0.0 if (np.isinf(dx) or dx <= 0) else 0.5 * 10 ** dose / dx
        else:                                             # CT99021, BIRB0796: no dose response in AGS
            h_cal = h_dir = max(0.0, fx) * 0.5; term_cal = term_dir = 0.0
        out["h" + key] = h_cal; out["hd" + key] = h_dir; out["t" + key] = term_cal; out["td" + key] = term_dir
    out["bliss_exp"] = 1 - (1 - out["hA"]) * (1 - out["hB"])
    out["bliss_exp_dir"] = 1 - (1 - out["hdA"]) * (1 - out["hdB"])
    out["bliss_excess"] = r.fAB - out["bliss_exp"]
    out["bliss_excess_dir"] = r.fAB - out["bliss_exp_dir"]
    out["hsa_margin_full"] = r.fAB - max(r.fA, r.fB)
    ok = r.fAB >= 0.1          # CI is undefined when the combination has (almost) no effect
    out["CI_cal"] = out["tA"] + out["tB"] if ok else np.nan; out["CI_dir"] = out["tdA"] + out["tdB"] if ok else np.nan
    return pd.Series(out)
res = pd.concat([df, df.apply(metrics, axis=1)], axis=1)
pd.set_option("display.width", 250)
cols = ["A", "B", "call", "control", "fA", "fB", "fAB", "hsa_margin_full", "bliss_exp", "bliss_excess", "bliss_excess_dir", "CI_cal", "CI_dir"]
print("\nPer combination (f = fractional growth inhibition 48 h after drug; control = cell-index increment of control):")
print(res[cols].round(2).to_string())
res.to_csv("out/flobak_synergy.tsv", sep="\t")
# classes
for m in ("hsa_margin_full", "bliss_excess", "bliss_excess_dir", "CI_cal", "CI_dir"):
    a = res.loc[res.call == "syn", m].dropna(); b = res.loc[res.call == "non", m].dropna()
    auc = stats.mannwhitneyu(a, b).statistic / (len(a) * len(b))
    print(f"{m:18} synergistic (n={len(a)}) median {a.median():+.2f} [{a.min():+.2f},{a.max():+.2f}] | non (n={len(b)}) median {b.median():+.2f} [{b.min():+.2f},{b.max():+.2f}] | P(syn>non) = {auc:.2f}")
# repeatability of the same single drug across experiments (full dose)
print("\nSame drug, same dose in different experiments: fractional inhibition (min..max, n)")
for d in ("PI103", "PD0325901", "5Z", "AKTi", "PKF118-310", "CT99021", "BIRB0796"):
    v = [r.fA for r in res.itertuples() if r.A == d] + [r.fB for r in res.itertuples() if r.B == d]
    print(f"  {d:11} {min(v):.2f}..{max(v):.2f} (n={len(v)}, sd {np.std(v, ddof=1):.2f})")

print("\nSensitivity of the Bliss excess to the half-dose/full-dose effect ratio r (same r for both drugs); the measured curves give r = 0.66 (PI103), 0.8 (PD0325901), 0.69 (5Z)")
for fig in (18, 21):
    r_ = res.loc[fig]
    print(f"  Fig {fig} {r_.A}+{r_.B}: " + ", ".join(f"r={q:.2f}: {r_.fAB - (1 - (1 - q*r_.fA)*(1 - q*r_.fB)):+.2f}" for q in (0.5, 0.66, 0.8, 0.9, 1.0)))
nonres = res.loc[res.call == "non", "bliss_excess"]
print(f"non-synergistic combinations: Bliss excess mean {nonres.mean():+.2f}, sd {nonres.std(ddof=1):.2f}; z of Fig 18 = {(res.loc[18,'bliss_excess']-nonres.mean())/nonres.std(ddof=1):.1f}, Fig 21 = {(res.loc[21,'bliss_excess']-nonres.mean())/nonres.std(ddof=1):.1f}")
