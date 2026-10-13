import pandas as pd, numpy as np, glob, os
for c in ['PI','PD','TAK','PITAK','PIPD','PIPD_vsPI','PIPD_vsPD','PIPD_dose','PIPD_add','PITAK_vsPI','PITAK_vsTAK','PITAK_dose']:
    row=[]
    for m in ['camera_fixed','camera_est','fry']:
        f=f'out/r_{m}_{c}.tsv'
        if not os.path.exists(f): row.append('NA'); continue
        d=pd.read_csv(f,sep='\t')
        row.append('%s p<.05 %d, BH<.05 %d'%(m,(d.PValue<.05).sum(),(d.q<.05).sum()))
    print('%-11s'%c,' | '.join(row))
