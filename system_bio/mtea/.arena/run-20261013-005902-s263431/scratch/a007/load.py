import pandas as pd, numpy as np, openpyxl
P='D:/Cander/School/3-2/system_bio/mtea/'
def moesm4():
    wb=openpyxl.load_workbook(P+'si/MOESM4.xlsx',data_only=False)
    out={}
    for ws in wb.worksheets:
        rows=list(ws.iter_rows(values_only=True))
        hdr=rows[0]
        data=[r for r in rows[2:] if r[0] is not None]
        df=pd.DataFrame([(int(r[0]),r[1],r[2],r[3]) for r in data],columns=['id','desc','score','p'])
        df['hdrcell']=hdr[2]
        out[ws.title]=df
    return out
