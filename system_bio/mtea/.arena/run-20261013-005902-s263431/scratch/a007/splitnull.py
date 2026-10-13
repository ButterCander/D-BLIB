import json, numpy as np, pandas as pd, itertools
from scipy import stats
from scipy.stats import false_discovery_control
P='D:/Cander/School/3-2/system_bio/mtea/'
rng=np.random.default_rng(11)
cnt=pd.read_csv(P+'GSE285616_counts.csv.gz',index_col='Symbol'); cnt=cnt[~cnt.index.duplicated()]; cnt=cnt.loc[cnt.sum(axis=1)>=10]
print(cnt.columns.tolist())
x=cnt.values.astype(float); pos=(x>0).all(axis=1)
sf=np.median(x[pos]/np.exp(np.log(x[pos]).mean(axis=1,keepdims=True)),axis=0)
Y=np.log2(x/sf+1); G=Y.shape[0]; genes=list(cnt.index); gi={g:i for i,g in enumerate(genes)}
batch=np.array([i//6 for i in range(24)]); cond=np.array([i%6 for i in range(24)])
tasks=json.load(open(P+'out/task_genes.json'))
sets={}
for t,gl in tasks.items():
    idx=sorted({gi[g] for g in gl if g in gi})
    if 5<=len(idx)<=200: sets[t]=np.array(idx)
names=list(sets); sizes=sorted({len(v) for v in sets.values()})
print('task sets',len(names),'genes',G)
def fit(Yd,X):
    XtXi=np.linalg.inv(X.T@X); beta=XtXi@X.T@Yd.T
    R=Yd.T-X@beta; df=24-np.linalg.matrix_rank(X); s2=(R**2).sum(0)/df
    return beta,XtXi,s2,df
def pvals_perm(z,ndraw=2000):
    n=len(z); out=np.empty(len(names))
    null={m:z[rng.integers(0,n,size=(ndraw,m))].mean(1) for m in sizes}
    for k,nm in enumerate(names):
        idx=sets[nm]; m=len(idx); zin=z[idx].mean(); nl=null[m]
        out[k]=(1+np.sum(np.abs(nl-nl.mean())>=abs(zin-nl.mean())))/(ndraw+1)
    return out
fpr=[];bh=[];rows=[]
for c in range(6):
    for split in [(0,1),(0,2),(0,3)]:
        a=[b for b in range(4) if b in split]; b_=[b for b in range(4) if b not in split]
        if 0 not in a: continue
        # new factor: condition c split into c_a (batches in a) and c_b
        lev=cond.copy()
        grp=np.zeros(24,int)
        # build design: intercept, batch(3), condition dummies for 5 others, c_a dummy, c_b dummy -> drop one level
        cols=[np.ones(24)]+[(batch==b).astype(float) for b in (1,2,3)]
        for k in range(6):
            if k==c: continue
            cols.append((cond==k).astype(float))
        isa=((cond==c)&np.isin(batch,a)).astype(float)
        cols.append(isa)
        X=np.column_stack(cols)
        beta,XtXi,s2,df=fit(Y,X)
        cvec=np.zeros(X.shape[1]); cvec[-1]=1
        est=cvec@beta; se=np.sqrt(s2*(cvec@XtXi@cvec)); t=est/se
        z=stats.norm.ppf(np.clip(stats.t.cdf(t,df),1e-12,1-1e-12))
        p=pvals_perm(z)
        fpr.append(np.mean(p<0.05)); bh.append(int((false_discovery_control(p,method='bh')<0.05).sum()))
        rows.append((c,split,df,fpr[-1],bh[-1]))
        print(c,split,'df',df,'FPR@.05=%.3f'%fpr[-1],'BH rej',bh[-1],'n sig genes p<.05 (z)',int((np.abs(z)>1.96).sum()),flush=True)
print('mean FPR',np.mean(fpr),'min',np.min(fpr),'max',np.max(fpr),'mean BH',np.mean(bh))
