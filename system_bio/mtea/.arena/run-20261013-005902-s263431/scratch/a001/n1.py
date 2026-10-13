import json, numpy as np, pandas as pd
g=json.load(open('out/eflux_grid.json'))
vals=[v for x in g for v in x['pred'].values()]
print('E-flux settings',len(g),'pred min %.2f max %.2f'%(min(vals),max(vals)))
rh=[x['rho'] for x in g]; print('rho range',min(rh),max(rh))
# which settings give growth zero?
for x in g[:3]: print({k:x[k] for k in ('vmax','floor','conc')},x['pred'])
# per-setting growth mean (is pred in percent? )
print([ (x['vmax'],x['floor'],round(np.mean(x['growth']),3)) for x in g])
f=pd.read_csv('out/flobak_synergy.tsv',sep='\t')
pd.set_option('display.width',250)
print(f[['fig','A','B','call','fA','fB','fAB','hsa_margin_full','bliss_excess','bliss_excess_dir','CI_cal','CI_dir']].round(3).to_string())
syn=f[f.call=='syn']; non=f[f.call!='syn']
print('n syn',len(syn),'non',len(non))
print('HSA syn range',syn.hsa_margin_full.min(),syn.hsa_margin_full.max(),'non',non.hsa_margin_full.min(),non.hsa_margin_full.max())
print('Bliss excess syn range',syn.bliss_excess.min(),syn.bliss_excess.max(),'non',non.bliss_excess.min(),non.bliss_excess.max())
print('non Bliss excess mean %.3f sd %.3f'%(non.bliss_excess.mean(),non.bliss_excess.std()))
print('CI syn',syn.CI_cal.min(),syn.CI_cal.max(),'non',non.CI_cal.dropna().min(),non.CI_cal.dropna().max(),'n non with CI',non.CI_cal.notna().sum())
from sklearn.metrics import roc_auc_score
