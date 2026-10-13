import pandas as pd, numpy as np, json
pd.set_option('display.width',250); pd.set_option('display.max_columns',40)
f=pd.read_csv('out/flobak_synergy.tsv',sep='\t',index_col=0)
print(f.columns.tolist()); print(f.round(3).to_string())
g=json.load(open('out/eflux_grid.json'))
print(type(g),len(g)); print(g[0] if isinstance(g,list) else list(g.items())[:1])
