from h import *
import json
x = pd.ExcelFile(P("si","MOESM5.xlsx"))
md = pd.read_excel(P("si","MOESM3.xlsx"), sheet_name="Task Metadata").set_index("id")
for sh in ["PI3Ki-TAKi","PI3Ki-MEKi"]:
    d=x.parse(sh).iloc[:,:6]; d.columns=["task","desc","A","B","AB","dir"]
    m = d[d.desc.str.contains("ydroxymethyl|cetone|cetoacet|HMG", case=False)]
    print(sh, m.round(2).values.tolist())
print(md[md.DESCRIPTION.str.contains("ydroxymethyl|HMG", case=False)][["SUBSYSTEM","DESCRIPTION","N.GENES"]])
tg = json.load(open(P("out","task_genes.json")))
for t in ["X138","X139","X140","X141","X142"]:
    print(t, md.loc[t,"DESCRIPTION"], tg.get(t))
