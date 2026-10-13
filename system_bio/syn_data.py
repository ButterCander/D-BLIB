"""Load LogicForests data (Taoma 2025). pickle is restricted to numpy so a malicious .pkl cannot run code."""
import pickle
import numpy as np
import pandas as pd

D = r"D:\Cander\School\3-2\system_bio\LogicForests\file"


class NumpyOnly(pickle.Unpickler):
    def find_class(self, module, name):
        if module.split(".")[0] == "numpy":
            return super().find_class(module.replace("numpy.core", "numpy._core"), name)
        raise pickle.UnpicklingError(f"blocked {module}.{name}")


def _load(f):
    return NumpyOnly(open(f"{D}\\{f}", "rb")).load()


X, y, proteins = _load("feat_mean.pkl"), _load("ave_hsa.pkl"), _load("name_list.pkl")
meta = pd.read_csv(f"{D}\\drugcomb_RF_574.csv", index_col=0)
