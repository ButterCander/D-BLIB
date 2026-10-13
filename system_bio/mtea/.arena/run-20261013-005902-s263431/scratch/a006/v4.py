import numpy as np
from PIL import Image
def grid_rows(im):
    # dotted horizontal gridlines: rows with many dark-gray pixels, not the frame
    g=(im.sum(axis=2)<420)&(abs(im[:,:,0]-im[:,:,1])<20)&(abs(im[:,:,1]-im[:,:,2])<20)
    cnt=g[:,60:850].sum(1)
    rows=[r for r in range(im.shape[0]) if cnt[r]>150]
    # merge adjacent
    out=[];
    for r in rows:
        if not out or r-out[-1][-1]>1: out.append([r])
        else: out[-1].append(r)
    return [np.mean(o) for o in out]
def colmask(im,name):
    R,G,B=im[...,0],im[...,1],im[...,2]
    return {'red':(R>200)&(G<70)&(B<70),'green':(R<90)&(G>190)&(B<90),'blue':(R<70)&(G<70)&(B>190),'magenta':(R>190)&(G<70)&(B>190)}[name]
def curve_y(im,name,xs):
    m=colmask(im,name)[:,xs[0]:xs[1]+1]
    ys=np.where(m.any(axis=1))[0]
    return ys
for fn,cols in [('si/fl_p18_img0.png',{'PD0325901':'green','PI103':'blue'}),('si/fl_p18_img3.png',{'PI103':'green','5Z':'blue'})]:
    im=np.array(Image.open(fn).convert('RGB')).astype(int)
    gr=grid_rows(im); print(fn,'grid rows',[round(x,1) for x in gr])
    # find x of end: x range last 8 px of plot
    for nm in ['red','green','blue','magenta']:
        ys=curve_y(im,nm,(836,848)); 
        print(nm, 'y range at t~69-70h:', ys.min() if len(ys) else None, ys.max() if len(ys) else None, 'median', np.median(ys) if len(ys) else None)
