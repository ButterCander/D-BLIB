import numpy as np, pandas as pd, pickle, itertools
from scipy import stats
res=pickle.load(open('.arena/run-20261013-005902-s263431/scratch/a002/w/ind_res.pkl','rb'))
SREBP=["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR=["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
def prog(l,genes=SREBP): return l.reindex([g for g in genes if g in l.index]).dropna().mean()
# per-line list of dataset-level program means (documented composite: mean over datasets)
D={'CAPAN2':[('tram24','CAPAN2')],'CFPAC1':[('tram24','CFPAC1'),('kras24','CFPAC1')],'PANC1':[('tram24','PANC1'),('kras24','PANC1')],
   'MIAPACA2':[('tram24','MiaPaca2'),('tram48','MiaPaca2')],'KP4':[('tram48','KP4')],'ASPC1':[('kras24','AsPC1')],'SW1990':[('kras24','SW1990')]}
ind={k:np.mean([prog(res[x]) for x in v]) for k,v in D.items()}
print('composite inducibility (recomputed):',{k:round(v,2) for k,v in ind.items()})
model=pd.read_csv('external/Model_24Q4.csv').set_index('ModelID')
n2i={r.StrippedCellLineName:i for i,r in model.iterrows() if r.OncotreeLineage=='Pancreas'}
E=pd.read_csv('external/expr_pancreas.tsv',sep='\t',index_col=0)
print('expr shape',E.shape, 'lines found',[k for k in D if k in n2i])
ind=pd.Series(ind)
def rho(pred): 
    x=ind[pred.index]; return stats.spearmanr(x,pred)
bt=pd.Series({k:E.loc[n2i[k],TR].mean() for k in D})
print('baseline trunk mean log2TPM:',bt.round(2).to_dict())
r=rho(bt); print('A) composite vs baseline trunk6: rho=%.3f p=%.4f'%(r[0],r[1]))
# exact permutation p
rk=stats.rankdata(ind[bt.index].values); rb=stats.rankdata(bt.values)
cnt=sum(abs(np.corrcoef(rk,np.array(p))[0,1])>=abs(np.corrcoef(rk,rb)[0,1])-1e-9 for p in itertools.permutations(rb)) ; print(' exact two-sided permutation p (Pearson on ranks) =',cnt/5040)
# leave-one-out
for k in D:
    s=[x for x in D if x!=k]; print(' LOO drop',k,'rho=%.2f'%stats.spearmanr(ind[s],bt[s])[0], end=';')
print()
# B) trametinib 24h only n=4
tr={'CAPAN2':prog(res[('tram24','CAPAN2')]),'CFPAC1':prog(res[('tram24','CFPAC1')]),'MIAPACA2':prog(res[('tram24','MiaPaca2')]),'PANC1':prog(res[('tram24','PANC1')])}
x=[bt[k] for k in tr]; y=list(tr.values()); print('B) trametinib24h only n=4 rho=%.2f'%stats.spearmanr(x,y)[0], {k:round(v,2) for k,v in tr.items()})
# C) outcomes: HMGCS1 only / trunk only (composite)
for nm,genes in (('HMGCS1 only',['HMGCS1']),('trunk6',TR),('program minus trunk',[g for g in SREBP if g not in TR])):
    i2=pd.Series({k:np.mean([prog(res[x],genes) for x in v]) for k,v in D.items()})
    print('C)',nm,'composite rho vs baseline trunk6 = %.2f'%stats.spearmanr(i2[bt.index],bt)[0],{k:round(v,2) for k,v in i2.items()})
# D) predictor not from trunk: baseline of other program genes / SREBF2 / INSIG1 / SCAP
other=[g for g in SREBP if g not in TR]
for nm,genes in (('program minus trunk',other),('INSIG1',['INSIG1']),('SREBF2',['SREBF2']),('SCAP',['SCAP']),('SQLE',['SQLE']),('HMGCS1 only',['HMGCS1']),('HMGCR only',['HMGCR'])):
    g=[x for x in genes if x in E.columns]
    b=pd.Series({k:E.loc[n2i[k],g].mean() for k in D}); print('D) baseline',nm,'rho=%.2f'%stats.spearmanr(ind[b.index],b)[0])
# E) distribution of |rho| for all genes as predictor (how special is 0.96 among 17k genes?)
Ei=E.loc[[n2i[k] for k in D]]; Ei.index=list(D)
allr=Ei.apply(lambda c: stats.spearmanr(ind[c.index],c)[0] if c.nunique()>3 else np.nan)
print('E) all genes: n=%d, |rho|>=0.96: %d (%.2f%%), >=0.9: %d (%.2f%%)'%(allr.notna().sum(),(allr.abs()>=.96).sum(),100*(allr.abs()>=.96).mean(),(allr.abs()>=.9).sum(),100*(allr.abs()>=.9).mean()))
# random 6-gene means
rng=np.random.default_rng(1); genes=allr.dropna().index.values
hits=0; N=3000
for _ in range(N):
    g=rng.choice(genes,6,replace=False); b=Ei[g].mean(axis=1); hits+= abs(stats.spearmanr(ind[b.index],b)[0])>=0.96-1e-9
print('   random 6-gene sets reaching |rho|>=0.96: %.2f%%'%(100*hits/N))
# among genes of other cholesterol-related? SREBP target genes' rank
print(' rho of 24 program genes baseline each:',{g:round(allr.get(g,np.nan),2) for g in SREBP if g in allr.index})
