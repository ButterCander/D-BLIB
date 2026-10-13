import pandas as pd, numpy as np
from scipy import stats
SREBP = ["HMGCS1","HMGCR","MVK","PMVK","MVD","IDI1","FDPS","FDFT1","SQLE","LSS","CYP51A1","MSMO1","NSDHL","HSD17B7","DHCR7","DHCR24","SC5D","EBP","TM7SF2","ACAT2","INSIG1","LDLR","SCD","FASN"]
TR = ["HMGCS1","HMGCR","MVK","MVD","IDI1","FDPS"]
res = {}
for ln in ['CAPAN2','CFPAC1','MiaPaca2','PANC1']:
    d = pd.read_csv(f'pdac/GSE238198_Pdac.{ln}.genes.normalizedCounts.txt.gz', sep='\t', index_col=0)
    d = d[~d.index.duplicated()]
    y = np.log2(d+1)
    nt = [c for c in y.columns if '_NT_' in c]
    out = {}
    for tp in ['3h','12h','24h']:
        tc = [c for c in y.columns if f'_{tp}_' in c]
        lfc = y[tc].mean(axis=1)-y[nt].mean(axis=1)
        out[tp]=lfc
        if tp=='24h':
            t = stats.ttest_ind(y[tc].values, y[nt].values, axis=1, equal_var=False)
            p = pd.Series(t.pvalue, index=y.index)
    L = pd.DataFrame(out)
    res[ln]=L
    prog = L.reindex([g for g in SREBP if g in L.index])
    print(ln, 'n genes', len(L), '| program mean (3/12/24h):', prog.mean().round(2).to_dict(), '| n program genes', len(prog), '| 24h up count', int((prog['24h']>0).sum()))
    print('   24h key genes:', {g: round(L.loc[g,'24h'],2) for g in ['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SQLE','LDLR','INSIG1']}, ' p(HMGCR)=%.3g'%p['HMGCR'])
    print('   mean counts NT:', {g: round(d[nt].loc[g].mean()) for g in ['HMGCS1','HMGCR','MVK','MVD','IDI1','SCD']})
# count lines with log2FC>=0.5 per gene at 24h
tab = pd.DataFrame({ln: res[ln]['24h'] for ln in res})
print(tab.loc[['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD','SREBF1','SQLE','LSS','DHCR7']].round(2))
print('lines with lfc>=0.5:', (tab.loc[['HMGCS1','HMGCR','MVK','MVD','IDI1','FDPS','SCD']]>=0.5).sum(axis=1).to_dict())
tab.to_pickle('.arena/run-20261013-005902-s263431/scratch/a004/pdac24.pkl')
