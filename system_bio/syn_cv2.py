"""Unseen-drug generalisation: leave-one-drug-out (train on samples without drug d, test on samples with d)."""
from syn_cv import *

tg = pd.read_csv("BoolSyn/data/drug_target.csv", encoding="utf-8-sig")
tg["drug"] = tg.drug.str.strip().replace({"ERLOTINIB": "ERLOTINIB HYDROCHLORIDE"})
print("drugs without target info:", [d for d in drugs if d not in set(tg.drug)])
targets = sorted(tg.protein.unique())
DT = np.array([[(tg.protein[tg.drug == d] == t).any() for t in targets] for d in drugs], float)  # drug x target
TARGET = np.minimum(DRUG @ DT, 1)  # sample x target (union of the two drugs' targets)

FEATS2 = {"Boolean(117)": X, "drug id+cell": np.hstack([DRUG, CELL]), "target+cell": np.hstack([TARGET, CELL]),
          "Boolean+target+cell": np.hstack([X, TARGET, CELL]), "Boolean+drug id+cell": np.hstack([X, DRUG, CELL])}
M2 = {k: MODELS[k] for k in ("Ridge", "SVR", "ExtraTrees", "HistGB")}


def drug_out(F, mk, min_n=15):
    r, a = [], []
    for i in range(len(drugs)):
        te = DRUG[:, i] > 0
        if te.sum() < min_n:
            continue
        p = mk().fit(F[~te], y[~te]).predict(F[te])
        r.append(pearsonr(y[te], p)[0] if p.std() > 0 else 0.0)
        a.append(roc_auc_score(y[te] > 0, p) if len(set(y[te] > 0)) == 2 else np.nan)
    return np.nanmean(r), np.nanmean(a), len(r)


if __name__ == "__main__":
    print(f"== leave-one-drug-out (mean over held-out drugs; Pearson r / AUC)")
    print(f"{'features':22}" + "".join(f"{m:>15}" for m in M2))
    for fname, F in FEATS2.items():
        row = ""
        for mk in M2.values():
            r, a, n = drug_out(F, mk)
            row += f"{r:8.3f}/{a:.3f}"
        print(f"{fname:22}{row}")
    print(f"(held-out drugs with >=15 samples: {n})")
