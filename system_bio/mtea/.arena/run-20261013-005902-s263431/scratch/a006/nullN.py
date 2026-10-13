import sys,time,numpy as np,pandas as pd
sys.path.insert(0,'mteapy/src'); sys.stdout.reconfigure(encoding='utf8')
from scipy.stats import false_discovery_control as bh
from cobra.io import read_sbml_model
from mteapy.tide import compute_TIDE
D="mteapy/src/mteapy/data"
gm=pd.read_csv(f"{D}/HumanGEM_genes.tsv",sep="\t").drop_duplicates("geneSymbol").set_index("geneSymbol")["geneID"]
ts=pd.read_csv(f"{D}/task_structure_matrix.tsv",sep="\t",index_col=0)
model=read_sbml_model(f"{D}/HumanGEM.xml.gz")
r0=pd.read_csv("out/de_PI103.tsv",sep="\t",index_col=0); r0["padj"]=r0["padj"].fillna(1)
mg=r0.loc[r0.index.intersection(gm.index)].copy(); mg.index=[gm[i] for i in mg.index]
mg=mg[~mg.index.duplicated()]
mg["lfc"]=np.where(mg.padj<.05,mg.lfc_shrunk,0.0)
au=pd.read_excel("si/MOESM4.xlsx",sheet_name="PI3Ki",header=None).iloc[2:,[0,2,3]]
au.columns=["tid","sa","pa"]; au["tid"]="X"+au.tid.astype(int).astype(str); au["pa"]=pd.to_numeric(au.pa); au["sa"]=pd.to_numeric(au.sa)
print("model genes with data",len(mg),"nonzero",(mg.lfc!=0).sum(),flush=True)
for N in [len(mg),3000,5000,16932]:
    e=pd.DataFrame({"lfc":mg.lfc.values},index=mg.index)
    extra=max(0,N-len(e))
    if extra: e=pd.concat([e,pd.DataFrame({"lfc":np.zeros(extra)},index=[f"dummy{i}" for i in range(extra)])])
    t0=time.time()
    res=compute_TIDE(e,"lfc",ts,model,"absmax",1000,1,True,7)
    nulls={r.task_id:np.array(r.random_score_array.split(";"),dtype=float) for r in res.itertuples()}
    p=[min((nulls[r.tid]<=r.sa).mean(),(nulls[r.tid]>=r.sa).mean()) if r.tid in nulls else np.nan for r in au.itertuples()]
    au["pu"]=p; ok=au.pu.notna()
    q=bh(np.clip(2*au.pu[ok].values,1e-12,1))
    print(f"N={N:6d} nonzero frac {(e.lfc!=0).mean():.2f}: authors' PI3Ki scores -> p<.025 {(au.pu[ok]<.025).sum()} of {ok.sum()}  BH<.05 {(q<.05).sum()}   [authors' own: 78 / 59]  {time.time()-t0:.0f}s",flush=True)
