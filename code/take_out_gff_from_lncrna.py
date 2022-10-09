import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from time import sleep
import requests
import json
from pathlib import Path
from tqdm import tqdm
from bs4 import BeautifulSoup
from tqdm import tqdm

def get_ncbi_code(nscaffold):
    URL = f"https://www.ncbi.nlm.nih.gov/nuccore/?term={nscaffold}"
    page = requests.get(URL)
    soup = BeautifulSoup(page.content, 'html.parser')
    if len(soup.select('.itemid')) > 0:
        results = soup.select('.itemid')[0]
        r = results.text[25:]
    else:
        r = "NOT_FOUND"
    return r

def df2dic(al, intervals):
    ranges = { interval : defaultdict(list) for interval in intervals }
    for Var in intervals:
        temp_df = al[al["type"]==Var]
        for index, row in temp_df.iterrows():
            ranges[Var][row['chr']].append([row['start'],row['stop']])
    return ranges

#quick check if two regions overlap
def check_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0])) > 0

def range_takeout(contig, evil):
    A0, A1 = contig #region that we like, e.g. conserved region, start and end respectively, start being the lower number no matter the direction of the gene for example
    B0, B1 = evil #region that we don't want in the end, e.g. coding region, start and end
    if( A1 < B0 or B1 < A0 ): # no overlap
        return [contig]
    elif ( A0 < B0 and B0 <= A1 and A1 <= B1) : # tail left
        return [[A0, B0]]
    elif (B0 <= A0 and A0 <= B1 and B1 < A1): # tail right
        return [[B1, A1]]
    elif ( B0 <= A0 and A1 <= B1): # superset
        return []
    elif (A0 < B0 and B1 < A1): # subset
        return [[A0,B0],[B1,A1]]
    else:
        print("ERROR:",contig, evil)

def subtract_scaffold(contig, enemy, debug=False, max_safety=1e6):
    good_idx = 0
    evil_idx = 0
    safety = 0 #this is to avoid an infinite while loop

    while True and safety < max_safety:
        if(good_idx >= len(contig) or evil_idx >= len(enemy)):
            print("YOU ARE A FAILURE.")
            break
        good = contig[good_idx]
        evil = enemy[evil_idx]
        good_start, good_end = good[0], good[1]
        evil_start, evil_end = evil[0], evil[1]
        if debug:
            print(f"checking {good} {evil} :",end=" ")
        if(check_overlap(good,evil)):
            #print("OVERLAPPED")
            #print(good, evil)
            if debug:
                print("overlap")
            segments = range_takeout(good,evil)
            if debug:
                print("segments",segments)
            contig = contig[:good_idx] + segments + contig[good_idx+1:]
            if debug:
                print("new",contig[good_idx-1:good_idx+2], "current idx", contig[good_idx])

            if(good_idx >= len(contig) or evil_idx >= len(enemy)):
                print("extra check saved us.")
                break

        else:
            if(evil_end <= good_start):
                if debug:
                    print("evil")
                evil_idx += 1
                if(evil_idx >= len(enemy)):
                    if debug:
                        print("EVIL FINISHED")
                    break
                if debug:
                    print("next evil", enemy[evil_idx])
            else:
                if debug:
                    print("good")
                good_idx += 1
                if(good_idx >= len(contig)):
                    if debug:
                        print("GOOD FINISHED")
                    break
                if debug:
                    print("next good", contig[good_idx])

        safety += 1
    if debug:
        print(safety)
    return contig

def get_sub_results(all_chrs, al, kind1, kind2, ranges):
    ref=ranges[kind1]
    alt=ranges[kind2]
    results = pd.DataFrame()
    for ch in all_chrs:
        if len(al[(al["chr"]==ch) & (al["type"]==kind1)]) > 0:
            if len(al[(al["chr"]==ch) & (al["type"]==kind2)]) > 0: # there is something to substract
              new_contig = subtract_scaffold(ref[ch], alt[ch])
              for c in new_contig:
                  results = results.append({'chr' : ch, 'start' : c[0], 'stop' : c[1]},
                          ignore_index = True)
            else:
              for index, row in al[(al["chr"]==ch) & (al["type"]==kind1)].iterrows():
                  results = results.append({'chr' : ch, 'start' : row.start, 'stop' : row.stop},
                          ignore_index = True)

    results.start = results.start.astype(int)
    results.stop = results.stop.astype(int)
    return results

df=pd.read_csv("sp4.lncRNAs_overlap_processed_5.0.csv", sep=",")
df=df.drop(columns=['#feat_name', 'mapped_int','recip','asm_unit','coverage','source_sub_start','source_sub_stop','source_int', 'source_length', 'mapped_length', 'source_strand', 'mapped_strand'])
print(len(df))
unmapped=len(df) - len(df.dropna(axis=0))
print("number of loci that didn't map anywhere in the 5.0 version")
print(unmapped) # number of loci that didn't map anywhere in the 5.0 version
df=df.dropna(axis=0)

ndf=df.groupby(['source_id','source_start']).count()
print("number of loci that mapped more than once")
print(len(ndf[ndf["mapped_id"]>=2])) # number of loci that mapped more than once

df=df.drop_duplicates(subset=['source_id','source_start']) # only kept first match

#spurscaff_to_NW={}
#spurscaffs=df.mapped_id.unique()
#for i in tqdm(range(len(spurscaffs))):
  #NW=get_ncbi_code(spurscaffs[i])
  #spurscaff_to_NW[spurscaffs[i]] = NW
#with open("spurscaff_to_NW.json", "w") as outfile:
  #json.dump(spurscaff_to_NW, outfile)

with open('spurscaff_to_NW.json', 'r') as f: # saved most already to a file because get_ncbi_code takes a while to run
  spurscaff_to_NW = json.load(f)

checklist=[]
for i in df.mapped_id.unique():
  if i not in spurscaff_to_NW:
    checklist.append(i)

if len(checklist) >= 1: #checking if I have the scaffold info for every chromosome I need
  for i in tqdm(range(len(checklist))):
    spurscaff_to_NW[checklist[i]]=get_ncbi_code(checklist[i])
    sleep(0.1)

df['mapped_id']=df['mapped_id'].map(spurscaff_to_NW)

# test if all remapping to scaffold has worked
print(len(df[df.isnull().any(axis=1)])) # should be 0

df=df.rename(columns={"mapped_id": "chr", "mapped_start": "start","mapped_stop": "stop"})

gff=pd.read_csv("nonoverlapping_nopromoter_gffannotations.csv")
gff=gff.dropna()

reference = "lnc"
subtractables = ["gff"]
intervals = [reference] + subtractables

#processing dataframes
df["type"] = "lnc"
gff["type"] = "gff"
gff=gff.sort_values(by=["chr",'start'])
df=df.sort_values(by=["chr",'start'])
al = pd.concat([gff, df])
all_chrs=al.chr.unique()
al["length"] = al["stop"] - al["start"]

#CHECK IF THERE ARE NEGATIVE NUMBERS
min(al["length"]) < 0

#make it into a dictionary
ranges = df2dic(al, intervals)

results = get_sub_results(all_chrs, al, "lnc", "gff", ranges)

results["length"]=results["stop"]-results["start"]
print(results.length.sum())

results.to_csv("lncRNA-gff.csv")
