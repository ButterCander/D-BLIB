import json, os, pandas as pd
from cobra.io import read_sbml_model
D="mteapy/src/mteapy/data"
ts=pd.read_csv(f"{D}/task_structure_matrix.tsv",sep="\t",index_col=0).astype(bool)
model=read_sbml_model(f"{D}/HumanGEM.xml.gz")
gs=pd.read_csv(f"{D}/HumanGEM_genes.tsv",sep="\t").drop_duplicates("geneID").set_index("geneID")["geneSymbol"]
rx2g={r.id:[gs.get(g.id) for g in r.genes if gs.get(g.id)] for r in model.reactions if r.id in ts.index}
out={t:sorted({g for r in ts.index[ts[t]] for g in rx2g.get(r,[])}) for t in ts.columns}
json.dump(out,open("out/task_genes.json","w"))
print("tasks",len(out),"median genes",sorted(len(v) for v in out.values())[len(out)//2])
