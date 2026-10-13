import pandas as pd, numpy as np, openpyxl
from scipy.stats import spearmanr, pearsonr
wb=openpyxl.load_workbook('si/MOESM4.xlsx',read_only=True,data_only=True)
S={}
for ws in wb.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    df=pd.DataFrame([r[:4] for r in rows[2:]],columns=['id','desc','score','p']).dropna(subset=['id'])
    df['task']='X'+df.id.astype(int).astype(str); df['score']=pd.to_numeric(df.score); df['p']=pd.to_numeric(df.p)
    S[ws.title]=df.set_index('task')
mine={k:pd.read_csv(f'out/tide_{f}.tsv',sep='\t').set_index('task_id') for k,f in [('PI3Ki','PI103'),('MEKi','PD0325901'),('TAKi','5Z'),('PI3Ki-TAKi','PI103+5Z'),('PI3Ki-MEKi','PI103+PD')]}
drop=['X124','X125','X126','X127','X128']
print('excluding 5 tryptophan tasks')
for s in ['TAKi','MEKi']:
    for k in ['MEKi','TAKi']:
        a=S[s].score; b=mine[k].score
        idx=[t for t in a.index if t in b.index and t not in drop]
        print(' sheet',s,'vs my',k, round(pearsonr(a[idx],b[idx])[0],3), 'spearman all', round(spearmanr(a.reindex(b.index).fillna(0),b.fillna(0))[0],3))
d2=openpyxl.load_workbook('si/MOESM2.xlsx',read_only=True)
for ws in d2.worksheets:
    rows=list(ws.iter_rows(values_only=True))
    df=pd.DataFrame([r[:5] for r in rows[1:]],columns=['cond','geneID','geneSymbol','lfc','padj']).dropna(subset=['geneSymbol'])
    df=df.drop_duplicates('geneSymbol').set_index('geneSymbol')
    out={}
    for k,f in [('PI3Ki','PI103'),('MEKi','PD0325901'),('TAKi','5Z'),('PI3Ki-TAKi','PI103+5Z'),('PI3Ki-MEKi','PI103+PD')]:
        de=pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0)
        j=df.join(de[['log2FoldChange','lfc_shrunk']],how='inner')
        out[k]=round(pearsonr(j.lfc,j.lfc_shrunk)[0],3)
    print('MOESM2 sheet',ws.title,'n DEG',len(df),'padj max',df.padj.astype(float).max(),'|lfc| min',df.lfc.abs().min(),out)
