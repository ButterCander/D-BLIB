import json, numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
cnt=pd.read_csv('GSE285616_counts.csv.gz',index_col='Symbol'); cnt=cnt[~cnt.index.duplicated()]
print(cnt.shape, list(cnt.columns)[:12])
cnt=cnt.loc[cnt.sum(axis=1)>=10]; x=cnt.values.astype(float); pos=(x>0).all(axis=1)
sf=np.median(x[pos]/np.exp(np.log(x[pos]).mean(axis=1,keepdims=True)),axis=0)
Y=np.log2(x/sf+1); genes=list(cnt.index); gi={g:i for i,g in enumerate(genes)}
batch=np.array([i//6 for i in range(24)]); cond=np.array([i%6 for i in range(24)])
X=np.column_stack([np.ones(24)]+[(batch==b).astype(float) for b in (1,2,3)]+[(cond==c).astype(float) for c in range(1,6)])
beta=np.linalg.lstsq(X,Y.T,rcond=None)[0]; R=Y.T-X@beta   # 24 x G
Rz=(R-R.mean(0)); Rz=(Rz/ (np.linalg.norm(Rz,axis=0,keepdims=True)+1e-12)).T  # G x 24
tasks=json.load(open('out/task_genes.json'))
rhos=[];ms=[]
for t,gl in tasks.items():
    idx=sorted({gi[g] for g in gl if g in gi})
    if 5<=len(idx)<=200:
        C=Rz[idx]@Rz[idx].T; m=len(idx); rho=(C.sum()-m)/(m*(m-1)); rhos.append(rho); ms.append(m)
print('n sets',len(rhos),'median mean pairwise residual corr',np.median(rhos),'median size',np.median(ms))
rng=np.random.default_rng(0)
# random gene sets of same sizes
rr=[]
for m in ms:
    idx=rng.choice(len(genes),m,replace=False); C=Rz[idx]@Rz[idx].T; rr.append((C.sum()-m)/(m*(m-1)))
print('random sets median',np.median(rr))
print('VIF at median', 1+(np.median(ms)-1)*np.median(rhos))
