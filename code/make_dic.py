import pandas as pd
import json

df=pd.read_csv("SPU_LOC_echino.txt", sep="\t", header=None, names=["ID","LOC","name","gen","test"]) # same as GenePageGeneralInfo_AllGenes.txt

# building dictionary with SPU
SPU_LOC={}
for index,row in df.iterrows():
  newlocs=[]
  newspus=[]
  for coli in ["ID","LOC","name","gen","test"]: #columns where there can be gene name
    if "LOC" in str(row[coli]): #if you find LOC, there could still be SPU
      stuff = row[coli].split("|")
      if len(stuff) > 1:
        locs=[i for i in stuff if "LOC" in i]
        for l in locs:
          newlocs.append(l.split(" ")[-1])
        spus= [i for i in stuff if "SPU" in i]
        for s in spus:
          newspus.append(s)
      else:
        newlocs.append(row[coli].split(" ")[-1])

    if "SPU" in str(row[coli]): #if you find SPU, there could still be LOC
      x = row[coli].split("|")
      locs=[i for i in x if "LOC" in i]
      if len(locs)>0:
        for l in locs:
          newlocs.append(l.split(" ")[-1])
      y=[i for i in x if "SPU" in i]
      for s in y:
        newspus.append(s)

  for l in set(newlocs):
    for s in set(newspus):
      if l in SPU_LOC.keys():
        SPU_LOC[l].append(s)
      else:
        SPU_LOC[l]=[s]

with open('SPU_LOC.json', 'w') as fp:
    json.dump(SPU_LOC, fp)
