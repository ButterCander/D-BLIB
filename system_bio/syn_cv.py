"""Drug-synergy regression: paper-style CV vs leak-free CV, Boolean features vs drug identity, several models."""
import warnings
from scipy.stats import pearsonr
from sklearn.ensemble import ExtraTreesRegressor, HistGradientBoostingRegressor, RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import mean_squared_error, roc_auc_score
from sklearn.model_selection import GroupKFold, LeaveOneGroupOut, StratifiedKFold, cross_val_predict
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from syn_data import *

warnings.filterwarnings("ignore")
norm = lambda s: s.str.replace("-", "")
assert (norm(meta.findx).values == norm(pd.Series(proteins)).values).all(), "meta/X row mismatch"

cell = meta.cell_line_name
pair = meta[["drug_row", "drug_col"]].apply(lambda r: "|".join(sorted(r)), axis=1)
drugs = sorted(set(meta.drug_row) | set(meta.drug_col))
DRUG = np.array([[d in (a, b) for d in drugs] for a, b in zip(meta.drug_row, meta.drug_col)], float)
CELL = pd.get_dummies(cell).values.astype(float)
print(f"{len(y)} samples, {pair.nunique()} unique drug pairs, {len(drugs)} drugs, "
      f"{(pair.map(pair.value_counts()) > 1).mean():.0%} of samples share their pair with another cell line\n")

FEATS = {"Boolean(117)": X, "drug+cell": np.hstack([DRUG, CELL]), "Boolean+drug+cell": np.hstack([X, DRUG, CELL])}
MODELS = {
    "RF": lambda: RandomForestRegressor(300, max_depth=4, min_samples_leaf=16, max_samples=.7, n_jobs=-1, random_state=0),
    "ExtraTrees": lambda: ExtraTreesRegressor(300, max_depth=6, min_samples_leaf=8, n_jobs=-1, random_state=0),
    "HistGB": lambda: HistGradientBoostingRegressor(max_depth=3, learning_rate=.05, max_iter=150, min_samples_leaf=20, random_state=0),
    "Ridge": lambda: make_pipeline(StandardScaler(), RidgeCV(alphas=np.logspace(0, 4, 13))),
    "SVR": lambda: make_pipeline(StandardScaler(), SVR(C=3, epsilon=1.0, gamma="scale")),
}
strat = cell + (y > 0).astype(str)
SPLITS = {  # name -> list of (train, test)
    "paper-style (random)": list(StratifiedKFold(5, shuffle=True, random_state=0).split(X, strat)),
    "pair-grouped": list(GroupKFold(5).split(X, y, pair)),
    "leave-cell-line-out": list(LeaveOneGroupOut().split(X, y, cell)),
}

if __name__ == "__main__":
    for sname, cv in SPLITS.items():
        print(f"== CV: {sname}   (Pearson r / AUC, out-of-fold pooled)")
        print(f"{'features':20}" + "".join(f"{m:>15}" for m in MODELS))
        for fname, F in FEATS.items():
            row = ""
            for mk in MODELS.values():
                p = cross_val_predict(mk(), F, y, cv=cv)
                row += f"{pearsonr(y, p)[0]:8.3f}/{roc_auc_score(y > 0, p):.3f}"
            print(f"{fname:20}{row}")
        print()
