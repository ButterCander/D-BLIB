import numpy as np, pandas as pd
from PIL import Image
S='si/'
def load(f): return np.array(Image.open(S+f).convert('RGB')).astype(int)
COL={"red":lambda p:(p[...,0]>200)&(p[...,1]<70)&(p[...,2]<70),"green":lambda p:(p[...,0]<90)&(p[...,1]>190)&(p[...,2]<90),"blue":lambda p:(p[...,0]<70)&(p[...,1]<70)&(p[...,2]>190),"magenta":lambda p:(p[...,0]>190)&(p[...,1]<70)&(p[...,2]>190)}
PXH=11.3;X0=53.0;T_DRUG=21.0
def read(fn,cols,T):
    im=load(fn); xr=lambda t:int(round(X0+PXH*t))
    def row_at(color,t,hw=2):
        m=COL[color](im)[:,xr(t)-hw:xr(t)+hw+1]; ys=np.where(m.any(axis=1))[0]; return float(np.median(ys)) if len(ys) else np.nan
    r21=np.nanmedian([row_at(c,T_DRUG-0.8) for c in ('red','green','blue','magenta')])
    scale=(460.0-r21)/1.5; val=lambda row:1.0+(r21-row)/scale
    out={'control':val(row_at('red',T)),'combo':val(row_at('magenta',T))}
    for d,c in cols.items(): out[d]=val(row_at(c,T))
    return out
F={18:('fl_p18_img0.png',{'PD0325901':'green','PI103':'blue'}),21:('fl_p18_img3.png',{'PI103':'green','5Z':'blue'})}
for T in [37,45,53,69]:
    print('--- t =',T,'h on axis ( =',T-21,'h after drug )')
    for fig,(fn,cols) in F.items():
        v=read(fn,cols,T); ctrl=v['control']-1
        res={d:round(1-(v[d]-1)/ctrl,2) for d in cols}; res['combo']=round(1-(v['combo']-1)/ctrl,2)
        print(' fig',fig,res,' control increment %.2f'%ctrl)
