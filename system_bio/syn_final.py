"""A: repeated (5x 5-fold) CV with nested, pair-grouped hyper-parameter tuning. Saves out-of-fold predictions to syn_oof.npz."""
import re
from joblib import Parallel, delayed
from sklearn.linear_model import Ridge
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold
from syn_cv2 import *  # data, DRUG, CELL, TARGET, FEATS2, pair, cell, strat, y, X

R = 5
FS = {"Boolean": X, "drug+cell": FEATS2["drug id+cell"], "target+cell": FEATS2["target+cell"],
      "Boolean+target+cell": FEATS2["Boolean+target+cell"], "Boolean+drug+cell": FEATS2["Boolean+drug id+cell"]}
GRIDS = {
    "Ridge": (lambda: make_pipeline(StandardScaler(), Ridge()), {"ridge__alpha": np.logspace(-1, 4, 6)}),
    "SVR": (lambda: make_pipeline(StandardScaler(), SVR(epsilon=1.0)), {"svr__C": [1, 3, 10, 30], "svr__gamma": ["scale", "auto"]}),
    "RF": (lambda: RandomForestRegressor(200, max_samples=.7, n_jobs=1, random_state=0), {"max_depth": [3, 5], "min_samples_leaf": [8, 16]}),
}
key = lambda *a: "__".join(re.sub(r"\W", "_", s) for s in map(str, a))
SCEN = {
    "random": lambda s: StratifiedKFold(5, shuffle=True, random_state=s).split(X, strat),  # paper-style
    "pair": lambda s: StratifiedGroupKFold(5, shuffle=True, random_state=s).split(X, strat, pair),  # new drug pairs
}


def fit_fold(scen, rep, fs, mn, tr, te):
    mk, grid = GRIDS[mn]
    gs = GridSearchCV(mk(), grid, cv=GroupKFold(3), scoring="neg_mean_squared_error")
    gs.fit(FS[fs][tr], y[tr], groups=pair.values[tr])
    return scen, rep, fs, mn, te, gs.predict(FS[fs][te])


if __name__ == "__main__":
    jobs = [(sc, r, fs, mn, tr, te) for sc, sp in SCEN.items() for r in range(R) for tr, te in sp(r) for fs in FS for mn in GRIDS]
    print(len(jobs), "outer fits")
    oof = {key(sc, fs, mn): np.zeros((R, len(y))) for sc in SCEN for fs in FS for mn in GRIDS}
    for sc, r, fs, mn, te, p in Parallel(n_jobs=-1)(delayed(fit_fold)(*j) for j in jobs):
        oof[key(sc, fs, mn)][r, te] = p
    np.savez("syn_oof.npz", **oof)
    print("saved syn_oof.npz")
