import pandas as pd
import json
import csv

with open('SPU_LOC.json', 'r') as f: # saved most already to a file because get_scaff_num takes a while to run
  SPU_LOC = json.load(f)
#SPU_LOC

with open('all_locs.txt', 'r') as fd:
    reader = csv.reader(fd)
    locs=[]
    for row in reader:
      locs.append(row)

new_list=[]
not_found=[]
for l in locs:
  if l[0] in SPU_LOC:
    ss=SPU_LOC[l[0]]
    for s in ss:
      new_list.append(s)
  else:
    not_found.append(l[0])

print(len(not_found)/len(locs)*100) # percent LOCS not found corresponding SPU for

textfile = open("all_spus.txt", "w")
for element in new_list:
    textfile.write(element + "\n")
textfile.close()
