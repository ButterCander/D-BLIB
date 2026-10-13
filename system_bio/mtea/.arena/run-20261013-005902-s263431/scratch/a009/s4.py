import openpyxl, pandas as pd, numpy as np, sys
sys.stdout.reconfigure(encoding='utf-8')
wb=openpyxl.load_workbook('si/MOESM5.xlsx',read_only=True,data_only=True)
for ws in wb.worksheets:
    rows=[r for r in ws.iter_rows(values_only=True) if any(c is not None for c in r)]
    print('==',ws.title,len(rows))
    if ws.title=='Summary':
        for r in rows: print(r)
    else:
        df=pd.DataFrame([r[:6] for r in rows[1:]],columns=['task','desc','A','B','AB','dir'])
        print(df.dir.value_counts().to_dict())
        df['strong']=(df.A==0)&(df.B==0)
        print('strong',df.strong.sum(),'weak',(~df.strong).sum())
        if ws.title=='PI3Ki-MEKi':
            d=df[~df.strong]
            ext=np.where(d.AB<0, np.minimum(d.A,d.B)-d.AB, d.AB-np.maximum(d.A,d.B))
            print('weak margin median',np.median(ext),'frac<0.1',(ext<0.1).mean(),'n',len(ext))
        if ws.title=='PI3Ki-TAKi':
            d=df[~df.strong]
            ext=np.where(d.AB<0, np.minimum(d.A,d.B)-d.AB, d.AB-np.maximum(d.A,d.B))
            print('weak margin median',np.median(ext),'frac<0.1',(ext<0.1).mean(),'n',len(ext))
        print(df[df.strong][['task','desc','AB']].to_string())
