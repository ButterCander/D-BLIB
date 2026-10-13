from h import *
print("cond  n p<.025  | BH(2p)<.05  BH(2p)<.10 | BH(p)<.05 | p==0 count | min p")
for sh in auth:
    A = a_scores(sh)
    p = A['p'].values
    p2 = np.minimum(2*p,1)
    q2 = bh(p2); q1=bh(p)
    print(f"{sh:12s} {len(A)} {int((p<0.025).sum()):3d} | {int((q2<0.05).sum()):3d} {int((q2<0.10).sum()):3d} | {int((q1<0.05).sum()):3d} | {int((p==0).sum())} | {p.min()}")
# direction
for sh in auth:
    A=a_scores(sh); s=A[A.p<0.025]
    print(sh, "up",(s.score>0).sum(),"down",(s.score<0).sum(), "frac down %.3f"%((s.score<0).mean()))
tot_up=sum((a_scores(sh).query('p<0.025').score>0).sum() for sh in auth); tot_dn=sum((a_scores(sh).query('p<0.025').score<0).sum() for sh in auth)
print("all five: up",tot_up,"down",tot_dn, tot_dn/(tot_up+tot_dn))
three=["PI3Ki","PI3Ki-TAKi","PI3Ki-MEKi"]
u=sum((a_scores(sh).query('p<0.025').score>0).sum() for sh in three); d=sum((a_scores(sh).query('p<0.025').score<0).sum() for sh in three)
print("three PI3K conds: up",u,"down",d,d/(u+d))
# X76 asn
for sh in auth:
    A=a_scores(sh); print(sh,"X76",A.loc['X76'].to_dict())
