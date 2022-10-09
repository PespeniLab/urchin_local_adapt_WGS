#This is the code Lapo and I made to find all potential enhancer regions known to urchins. We do this by combining ATAC-seq, DNA-seq, Chip-seq, eRNA, and data about where L.var and S.pur are conserved.

#If an ATAC and Chip-seq region overlaps, we make it into one big region, and we continue to get all the ranges from all the rources. Once we combined everything mentioned above, we take out regions that overlap with coding regions (NCBI) and lncRNA (https://pubmed.ncbi.nlm.nih.gov/25959816/), in another code.

# ATAC-seq (all 4), DNAseq (all 2) and Chip-seq dfs didn't have any overlap within the df, however, 
# L.var did -> code to collapse those into single regions -> this code takes updated version.
# same with arenas_mena_lit and eRNA 

import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from time import sleep
import requests
import json
from pathlib import Path
from tqdm import tqdm

def get_segs(df, chromosome): #extracts pairs of start and stop for a chromosome from a df
    return [(start,stop) for _,(_,start,stop) in df[df['chr']==chromosome].iterrows()]

def check_overlap(a, b): #quick check if two regions overlap
    return max(0, min(a[1], b[1]) - max(a[0], b[0])) > 0

def clean_segments(src): #takes list of start stop, checks one against others, returns new (most of the time) extended region for overlaps
    if len(src) <= 1:
        return src
    a = src[0]
    dest = []
    for b in src[1:]:
        #print(a,b, check_overlap(a,b))
        if check_overlap(a,b):
            a = (min(a[0],b[0]),max(a[1],b[1]))
        else:
            dest.append(a)
            a = b
    dest.append(a)
    return dest

def merge_dfs(al, all_chrs):
  results = pd.DataFrame()
  i=1
  for ch in tqdm(range(len(all_chrs))):
      #print(i)
      src = get_segs(al,all_chrs[ch])
      tmp = clean_segments(src)
      for t in tmp:
          results = results.append({'chr' : all_chrs[ch], 'start' : t[0], 'stop' : t[1]}, ignore_index = True) 
      i=i+1
  return results

def process_results(df1, df2):
  #processing dataframes
  df1=df1.sort_values(by=["chr",'start'])
  df2=df2.sort_values(by=["chr",'start'])
  al = pd.concat([df1, df2])
  al["length"] = al["stop"] - al["start"]

  #CHECK IF THERE ARE NEGATIVE NUMBERS
  print(min(al["length"]) <= 0)

  al=al.drop(columns=['length'])
  al=al.sort_values(by=["chr",'start'])
  all_chrs=al.chr.unique()

  results=merge_dfs(al, all_chrs)

  return(results)

"""Overlap between different ATAC sources"""

atac1 = pd.read_csv("ATAC_all_24h.csv", skiprows=1)
atac1 = atac1.drop(columns=["Peak Name", "F-seq Peak Score"])
atac1 = atac1.rename(columns={"Scaffold": "chr", "Peak Start":"start", "Peak End": "stop"})

atac2 = pd.read_csv("ATAC_both_PMC_other_24h.csv", skiprows=1)
atac2 = atac2.drop(columns=["Peak Name"])
atac2 = atac2.rename(columns={"Scaffold": "chr", "Peak Start":"start", "Peak End": "stop"})

r1=process_results(atac1, atac2)
r1.to_csv("r1.csv")

atac3 = pd.read_csv("ATAC_other_24h.csv", skiprows=1)
atac3 = atac3.rename(columns={"Scaffold": "chr", "Peak Start":"start", "Peak End": "stop"})

r2=process_results(r1, atac3)
r2.to_csv("r2.csv")

atac4 = pd.read_csv("ATAC_PMC_24h.csv", skiprows=1)
atac4 = atac4.rename(columns={"Scaffold": "chr", "Peak Start":"start", "Peak End": "stop"})

r3=process_results(r2, atac4)
r3.to_csv("r3.csv")

r3=pd.read_csv("r3.csv")
print(r3.head())
r3=r3[["chr","start","stop"]]
r3["start"]= r3["start"].astype(int)
r3["stop"]= r3["stop"].astype(int)

"""Overlap between ATAC and Lvar"""

lvar=pd.read_csv("lytVar22_strPur31_48_50_overlap_processed.csv")
#lvar=lvar.rename(columns={"second.seqnames": "chr", "second.start": "start", "second.end": "stop"})
lvar=lvar[["chr","start","stop"]]

r4=process_results(r3, lvar)
r4.to_csv("r4.csv")

"""Overlap between that and Chip-seq"""

chip = pd.read_csv("Chip_all_peaks.csv", skiprows=31)
chip = chip[["chr","start","end"]]
chip=chip.rename(columns={"end": "stop"})

r5=process_results(r4, chip)
r5.to_csv("r5.csv")

"""Overlap between that and DNA-seqs"""

dna1=pd.read_csv("DNAseq_PMCminus.csv", skiprows=1)
dna1=dna1.rename(columns={"Scaffold":"chr","Start": "start", "End": "stop"})
dna1=dna1[["chr","start","stop"]]

r6=process_results(r5, dna1)
r6.to_csv("r6.csv")

dna2=pd.read_csv("DNAseq_PMCplus.csv", skiprows=1)
dna2=dna2.rename(columns={"Scaffold":"chr","Start": "start", "End": "stop"})
dna2=dna2[["chr","start","stop"]]

r7=process_results(r6, dna2)
r7.to_csv("r7.csv")

"""Overlap between that and arenas mena lit review"""

aren=pd.read_csv("arenas_mena_lit_overlap_processed.csv")
#aren=aren.rename(columns={"Sccaffold #": "chr", "end": "stop"})
aren=aren[["chr","start","stop"]]

r8=process_results(r7, aren)
r8.to_csv("r8.csv")
"""Overlap between that and eRNA"""

eRNA=pd.read_csv("eRNA_all_overlap_processed.csv")
#eRNA=eRNA.rename(columns={"Scaffold (genome v3.1)": "chr", "Start (genome v3.1)": "start", "End (genome v3.1)": "stop"})
eRNA=eRNA[["chr","start","stop"]]

r9=process_results(r8, eRNA)
r9.to_csv("overlapped_enhancers.csv")
print("DONE")

