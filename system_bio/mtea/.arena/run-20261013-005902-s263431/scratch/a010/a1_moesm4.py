import pandas as pd, numpy as np
from scipy import stats
base = r"D:\Cander\School\3-2\system_bio\mtea"
xl = pd.ExcelFile(base+r"\si\MOESM4.xlsx")
auth = {}
for sh in xl.sheet_names:
    raw = xl.parse(sh, header=None)
    print(sh, "header C1 =", raw.iloc[0,2], "| rows", raw.shape)
    df = xl.parse(sh, header=1)
    df.columns = [str(c) for c in df.columns]
    auth[sh] = df
    print(df.columns.tolist())
    n = len(df)
    p = df['p-value'].astype(float)
    sig = (p<0.025).sum()
    flag = df['Is Significant'].astype(str).str.lower().eq('true').sum()
    up = ((p<0.025)&(df['TIDE score']>0)).sum(); dn = ((p<0.025)&(df['TIDE score']<0)).sum()
    print("   n tasks", n, "p<0.025:", sig, "flag True:", flag, "up", up, "down", dn, "p<=0.025", (p<=0.025).sum())
import pickle
pickle.dump(auth, open("auth_moesm4.pkl","wb"))
