import pandas as pd
from bs4 import BeautifulSoup
import requests
from time import sleep
import numpy as np
from tqdm import tqdm
import json

df=pd.read_csv("input_for_extra_regannot_notannot_3.1.csv", sep=",")
df=df.drop(columns=['#feat_name', 'mapped_int','recip','asm_unit','coverage','mapped_stop','source_sub_start','source_sub_stop','source_int', 'source_length', 'mapped_length', 'source_strand', 'mapped_strand', 'source_stop'])

unmapped=len(df) - len(df.dropna(axis=0))
print("number of loci that didn't map anywhere in the 3.1 version")
print(unmapped) # number of loci that didn't map anywhere in the 3.1 version
df=df.dropna(axis=0)

ndf=df.groupby(['source_id','source_start']).count()
print("number of loci that mapped more than once")
print(len(ndf[ndf["mapped_id"]>=2])) # number of loci that mapped more than once

df=df.drop_duplicates(subset=['source_id','source_start']) # only kept first match


#NW to Scaffold number
def get_scaff_num(nscaffold):
  link=f"https://www.ncbi.nlm.nih.gov/nuccore/{nscaffold}?report=genbank"
  response=requests.get(link)
  soup=BeautifulSoup(response.text, 'html.parser')
  mystr=str(soup)
  splitted=mystr.split('Spur_3.1 ')
  if len(splitted) > 1:
    split2=splitted[2].split(',')
    return(split2[0])
  else:
    print("not found")
    print(nscaffold)
    return("NOT FOUND")

#with open("NW_to_scaf.json", "w") as outfile:
    #json.dump(NW_to_scaf, outfile)

with open('../supp_materials/NW_to_scaf_longest.json', 'r') as f: # saved most already to a file because get_scaff_num takes a while to run
  NW_to_scaf = json.load(f)

checklist=[]
for i in df.mapped_id.unique():
  if i not in NW_to_scaf:
    checklist.append(i)

if len(checklist) >= 1: #checking if I have the scaffold info for every chromosome I need
  for i in tqdm(range(len(checklist))):
    NW_to_scaf[checklist[i]]=get_scaff_num(checklist[i])
    sleep(0.1)

df['mapped_id']=df['mapped_id'].map(NW_to_scaf)

# test if all remapping to scaffold has worked
print(len(df[df.isnull().any(axis=1)])) # should be 0

"""Import stuff"""

#checked for each
#atac1['diff']=atac1['Peak End'] - atac1['Peak Start']
#if len(atac1[atac1["diff"]<=0]) > 0:
  #print("jajj")

ex_data=["arenas_mena_lit.csv", "Chip_all_peaks.csv", "ATAC_DNA_overlap.csv", "eRNA_all.csv", "lytVar22_strPur31_48_50.tsv", "lytVar22_strPur31_49_50.tsv", "lytVar22_strPur31_50_50.tsv", "sp4.lncRNAs.bed"]
tpath="../supp_materials/"
external_data=[tpath+i for i in ex_data]
print(external_data)

def check_lnc(mydf, theirdf, enhdf, label):
  for index, row in mydf.iterrows():
    #print(index)
    ch = row.mapped_id
    pos = row.mapped_start
    temp=theirdf[(theirdf['Scaffold']==ch) & (theirdf['Peak Start']<=pos) & (theirdf['Peak End']>=pos)]
    if len(temp) > 0:
      enhdf = enhdf.append({'source_id' : row.source_id, 'source_start' : row.source_start, 'mapped_id' : ch, 'mapped_start' : pos, 'region' : label },ignore_index = True)
    else:
      enhdf = enhdf.append({'source_id' : row.source_id, 'source_start' : row.source_start, 'mapped_id' : ch, 'mapped_start' : pos, 'region' : "not_lnc" },ignore_index = True)

  return enhdf

def check_enhancers(mydf, theirdf, enhdf, label):
  for index, row in mydf.iterrows():
    #print(index)
    ch = row.mapped_id
    pos = row.mapped_start
    temp=theirdf[(theirdf['Scaffold']==ch) & (theirdf['Peak Start']<=pos) & (theirdf['Peak End']>=pos)]
    if len(temp) > 0:
        enhdf = enhdf.append({'chr50' : row.source_id, 'pos50' : row.source_start, 'chr31' : ch, 'pos31' : pos, 'region' : label },ignore_index = True)
  return enhdf

def check_enhancers_2(mydf, theirdf, enhdf, label):
  for index, row in mydf.iterrows():
    #print(index)
    ch = row.source_id
    pos = row.source_start
    temp=theirdf[(theirdf['Scaffold']==ch) & (theirdf['Peak Start']<=pos) & (theirdf['Peak End']>=pos)]
    if len(temp) > 0:
        enhdf = enhdf.append({'chr50' : row.source_id, 'pos50' : row.source_start, 'chr31' : row.mapped_id, 'pos31' : row.mapped_start, 'region' : label },ignore_index = True)
  return enhdf

#lncRNA
outdf=pd.DataFrame()
for data in external_data:
  if "sp4.lncRNA" in data:
    print(data)
    names=["Scaffold", "Peak Start", "Peak End", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9"]
    temp=pd.read_csv(data, skiprows=0, sep="\t", names=names)
    outdf=check_lnc(df,temp, outdf, data)
    print("done")

lncdf=outdf[outdf["region"]!="not_lnc"] #separate lncRNA from other type of data, if a pos it in both ATAC and lncRNA, just keep it as lncRNA
lncdf.to_csv("lncdf.csv")
df=outdf[outdf["region"]=="not_lnc"]

#ATAC
enhdf=pd.DataFrame()
for data in external_data:
  if "ATAC" in data: # for now I won't differentiate between these
    print(data)
    temp=pd.read_csv(data, skiprows=1)
    enhdf=check_enhancers(df,temp, enhdf, "ATAC")
    print("done")

enhdf=enhdf.drop_duplicates(subset=['chr50', 'pos50']) # only keep 1 row for ATAC match

enhdf2=pd.DataFrame()
for data in external_data:
  if "lytVar22" in data:
    print(data)
    temp=pd.read_csv(data, skiprows=0, sep="\t")
    temp=temp.rename(columns={"second.seqnames": "Scaffold", "second.start": "Peak Start", "second.end": "Peak End"})
    enhdf2=check_enhancers(df,temp, enhdf2, data)
    print("done")
enhdf2=enhdf2.drop_duplicates(subset=['chr50', 'pos50'], keep='last') # only keep 1 row for lvar match
enhdf=pd.concat([enhdf, enhdf2]).reset_index()
enhdf=enhdf.drop(columns=['index'])

#other
for data in external_data:
  if "DNAseq" in data:
    print(data)
    temp=pd.read_csv(data, skiprows=1)
    temp=temp.rename(columns={"Start": "Peak Start", "End": "Peak End"})
    enhdf=check_enhancers(df,temp, enhdf, data)
    print("done")
  if "eRNA_all" in data:
    print(data)
    temp=pd.read_csv(data, skiprows=1)
    temp=temp.rename(columns={"Scaffold (genome v5.0)": "Scaffold", "Start (genome v5.0)": "Peak Start", "End (genome v5.0)": "Peak End"})
    temp.Scaffold = temp.Scaffold.str[:-2]
    enhdf=check_enhancers_2(df,temp, enhdf, data) # USEING ORIGINAL 5.0 version with this one
    print("done")
  if "arenas" in data:
    print(data)
    temp=pd.read_csv(data, skiprows=1)
    temp=temp.rename(columns={"Sccaffold #": "Scaffold", "start": "Peak Start", "end": "Peak End"})
    enhdf=check_enhancers(df,temp, enhdf, data)
    print("done")
  if "Chip_all" in data:
    print(data)
    temp=pd.read_csv(data, skiprows=31)
    temp=temp.rename(columns={"chr": "Scaffold", "start": "Peak Start", "end": "Peak End"})
    enhdf=check_enhancers(df,temp, enhdf, data)
    print("done")

#enhdf[enhdf.duplicated(subset=['chr50', 'pos50'], keep=False)].sort_values(by=["chr50","pos50"])
newdf=enhdf.groupby(["chr50","pos50", "chr31", "pos31"]).agg(lambda x: list(x)).reset_index()
#newdf['region'] = newdf['region'].apply(lambda x: x if x != ["ATAC", "sp4.lncRNAs.bed"] else ["sp4.lncRNAs.bed"])
newdf["confidence"]=newdf.region.str.len()
newdf.to_csv("reg_regions_withconfidence.csv")

