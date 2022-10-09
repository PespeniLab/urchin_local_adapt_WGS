import pandas as pd
import re

df=pd.read_csv("gathered_annotation.csv")
df=df[(df["Not_annot"]==0) & (df["pseudogene"]==0) & (df["snRNA"]==0) & (df["tRNA"]==0) & (df["lnc_RNA"]==0) & (df["snoRNA"]==0) & (df["rRNA"]==0) & (df["miRNA"]==0)]
df=df.drop(columns=['Not_annot', 'pseudogene','snRNA','tRNA','lnc_RNA','snoRNA','rRNA','miRNA','num_gene_hit','check','accunted_for'])

non_proms=df[df["promoter"]==0]
non_prom_locs=[]
for index, row in non_proms.iterrows():
  ls=row.gene_x.split(",")
  if len(ls) > 1:
    for l in ls:
      f=re.sub("[^0-9a-zA-Z]","", l)
      non_prom_locs.append(f)
  else:
    f=re.sub("[^0-9a-zA-Z]","", ls[0])
    #print(f)
    non_prom_locs.append(f)

proms=df[df["promoter"]!=0]
prom_overl_locs=[]
for index, row in proms.iterrows():
  ls=row.gene_x.split(",")
  if len(ls) > 1:
    for l in ls:
      f=re.sub("[^0-9a-zA-Z]","", l)
      prom_overl_locs.append(f)
  else:
    f=re.sub("[^0-9a-zA-Z]","", ls[0])
    #print(f)
    prom_overl_locs.append(f)

for p in prom_overl_locs:
  non_prom_locs.append(p)

proms=proms[proms["prom_region"]=="['protein_coding']"]
prom_locs=[]
for index, row in proms.iterrows():
  ls=row.prom_gene.split(",")
  if len(ls) > 1:
    for l in ls:
      f=re.sub("[^0-9a-zA-Z]","", l)
      prom_locs.append(f)
  else:
    f=re.sub("[^0-9a-zA-Z]","", ls[0])
    #print(f)
    prom_locs.append(f)

prom_locs=set(prom_locs)
non_prom_locs=set(non_prom_locs)
all_locs=set.union(prom_locs, non_prom_locs)

textfile = open("prom_locs.txt", "w")
for element in prom_locs:
    textfile.write(element + "\n")
textfile.close()

textfile = open("non_prom_locs.txt", "w")
for element in non_prom_locs:
    textfile.write(element + "\n")
textfile.close()

textfile = open("all_locs.txt", "w")
for element in all_locs:
    textfile.write(element + "\n")
textfile.close()
