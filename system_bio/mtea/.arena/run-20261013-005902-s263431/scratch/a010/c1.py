import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')
de={k:pd.read_csv(f'out/de_{f}.tsv',sep='\t',index_col=0) for k,f in {'5Z':'5Z','PD':'PD0325901','PI':'PI103','PI5Z':'PI103+5Z','PIPD':'PI103+PD'}.items()}
for s in ['TAKi','MEKi','PI3Ki','PI3Ki-TAKi','PI3Ki-MEKi']:
    d=pd.read_excel('si/MOESM2.xlsx',sheet_name=s).dropna(subset=['geneSymbol']).drop_duplicates('geneSymbol').set_index('geneSymbol')
    row=[s,len(d)]
    for k,x in de.items():
        j=d.join(x,how='inner',rsuffix='_m')
        row.append((k,len(j),round(np.corrcoef(j.log2FoldChange,j['log2FoldChange_m'])[0,1],3), round(np.corrcoef(j.log2FoldChange,j['lfc_shrunk'])[0,1],3)))
    print(row)
