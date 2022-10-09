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
                #print("extra check saved us.")
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
    for item in tqdm(range(len(all_chrs))):
        ch=all_chrs[item]
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

enhdf=pd.read_csv("overlapped_enhancers.csv")

lnc=pd.read_csv("sp4.lncRNAs_overlap_processed.csv")

"""Take out lnc"""

reference = "enh"
subtractables = ["lnc"]
intervals = [reference] + subtractables

#processing dataframes
enhdf["type"] = "enh"
lnc["type"] = "lnc"
enhdf=enhdf.sort_values(by=["chr",'start'])
lnc=lnc.sort_values(by=["chr",'start'])
al = pd.concat([lnc, enhdf])
all_chrs=al.chr.unique()
al["length"] = al["stop"] - al["start"]

#CHECK IF THERE ARE NEGATIVE NUMBERS
min(al["length"]) < 0

#make it into a dictionary
ranges = df2dic(al, intervals)

results = get_sub_results(all_chrs, al, "enh", "lnc", ranges)

results["length"]=results["stop"]-results["start"]
print(results.length.sum())

results.to_csv("enhancers-lncRNA.csv")

