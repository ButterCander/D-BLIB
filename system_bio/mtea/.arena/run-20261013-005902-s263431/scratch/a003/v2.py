import pandas as pd, numpy as np
x=pd.ExcelFile('si/MOESM5.xlsx'); print(x.sheet_names)
for s in x.sheet_names:
    d=pd.read_excel(x,s,header=None); print(s,d.shape); print(d.head(5).iloc[:,:8].to_string())
