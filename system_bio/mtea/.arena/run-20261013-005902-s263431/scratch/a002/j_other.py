import numpy as np, pandas as pd
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
def prog(d):
    x=d.reindex(SREBP).dropna(); return x.mean(), int((x>0).sum()), len(x)
# GSE206267
c=pd.read_csv("pdac/GSE206267_KP4_MiaPaca2_DMSO_Tram_rawcounts.txt.gz",sep="\t",index_col=0); c=c[c.sum(axis=1)>=30]
gi=pd.read_csv("tsirvouli/gene_info.tsv",sep="\t",index_col=0)["gene_name"]; c["sym"]=gi.reindex(c.index).values
c=c.dropna(subset=["sym"]).drop_duplicates("sym").set_index("sym"); nc=c/c.sum()*c.sum().median(); L=np.log2(nc+1)
for line,pre in (("KP4","KP4"),("MiaPaCa2","Mia")):
    d=L[[f"{pre}-T100-{i}" for i in (1,2,3)]].mean(axis=1)-L[[f"{pre}-DMSO-{i}" for i in (1,2,3)]].mean(axis=1)
    print("GSE206267",line,"prog %.2f up %d/%d"%prog(d), {g:round(d.get(g,np.nan),2) for g in ["HMGCS1","HMGCR","MVK","MVD"]})
# GSE279749
d=pd.read_csv("pdac/GSE279749_gene_fpkm_cell_lines.txt.gz",sep="\t").drop_duplicates("gene_name").set_index("gene_name")
for k,n in (("A","AsPC-1"),("P","PANC-1"),("C","CFPAC-1"),("S9","SW1990")):
    lf=np.log2(d[[f"{k}_M_{i}" for i in (1,2,3)]]+1).mean(axis=1)-np.log2(d[[f"{k}_NT_{i}" for i in (1,2,3)]]+1).mean(axis=1)
    print("GSE279749",n,"prog %.2f up %d/%d"%prog(lf), {g:round(lf.get(g,np.nan),2) for g in ["HMGCS1","HMGCR","MVK","MVD"]})
# GSE315130
for f,n,cols in (("PSN1","PSN1",("PSN1_6236","PSN1_Con")),("SW1990","SW1990",("SW1990_6236","SW1990_Con"))):
    try:
        t=pd.read_csv(f"pdac/GSE315130_{f}_cell_line_gene_sample_FPKM_with_symbol.txt.gz",sep="\t")
    except Exception as e:
        print(e); continue
    print(t.columns.tolist()[:12])
    t=t.drop_duplicates("Symbol").set_index("Symbol")
    tc=[c for c in t.columns if c.startswith(cols[0]+"_")]; cc=[c for c in t.columns if c.startswith(cols[1]+"_")]
    lf=np.log2(t[tc]+1).mean(axis=1)-np.log2(t[cc]+1).mean(axis=1)
    print("GSE315130",n,tc,cc,"prog %.2f up %d/%d"%prog(lf))
