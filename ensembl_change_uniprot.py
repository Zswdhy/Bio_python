import mygene

mg = mygene.MyGeneInfo()
xli = [
    "ENSG00000214562",
    "ENSG00000145113",
]
# scopes原始ID类型，fields转换后ID类型
out = mg.querymany(xli, scopes="ensembl.gene", fields="uniprot", species="human,mouse")
for item in out:
    print(item)
    # print("uniprotID", item["uniprot"]["Swiss-Prot"])
