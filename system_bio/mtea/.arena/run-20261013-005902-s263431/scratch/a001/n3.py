import pandas as pd, numpy as np
gs=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0,nrows=3)
print(gs.columns[:6].tolist(), gs.shape)
gs=pd.read_csv('external/crispr_gastric_rows.tsv',sep='\t',index_col=0)
cols=[c for c in gs.columns if c.startswith('HMGCS1') or c.startswith('HMGCR') or c.startswith('SQLE') or c.startswith('MVK') or c.startswith('GGPS1') or c.startswith('FNTA') or c.startswith('LSS') or c.startswith('ASNS')]
print(cols)
print('AGS ACH-000880' , 'ACH-000880' in gs.index)
print(gs.loc['ACH-000880',cols].round(2).to_dict())
print(gs[cols].median().round(2).to_dict())
