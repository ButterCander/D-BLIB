import pandas as pd, numpy as np
B="D:/Cander/School/3-2/system_bio/mtea/"
x = pd.ExcelFile(B+"si/MOESM5.xlsx")
print(x.sheet_names)
print(x.parse("Summary", header=None).to_string())
for sh in ["PI3Ki-TAKi","PI3Ki-MEKi"]:
    raw = x.parse(sh, header=None)
    print(sh, raw.shape)
    print(raw.iloc[:3,:8].to_string())
