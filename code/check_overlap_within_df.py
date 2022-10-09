# takes a supplementary regulatory region file and returns a dataframe for each region with overlaps resolved as one longer region
# so that overlapping positions and not counted multiple times
# run this code for eRNA, lncRNA, L.var, arenas_mena_lit review
# a separate code was run to process the gff file

import pandas as pd
from tqdm import tqdm

names=["chr", "start", "stop", "n1", "n2", "n3", "n4", "n5", "n6", "n7", "n8", "n9"]
lncRNA=pd.read_csv("sp4.lncRNAs.bed", skiprows=0, sep="\t", names=names)
atac1=lncRNA[["chr","start","stop"]]

atac1=atac1.sort_values(by=["chr",'start'])
all_chrs=atac1.chr.unique()

def clean_segments(src):
    if len(src) <= 1:
        return src

    a = src[0]
    dest = []
    for b in src[1:]:
        # print(a,b, check_overlap(a,b))
        if check_overlap(a,b):
            #print("OVERLAP")
            #print(src)
            a = (min(a[0],b[0]),max(a[1],b[1]))
        else:
            dest.append(a)
            a = b
    dest.append(a)
    return dest

#returns TRUE if there is overlap
def check_overlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0])) > 0
def get_segs(df, chromosome):
    return [(start,stop) for _,(_,start,stop) in df[df['chr']==chromosome].iterrows()]

results = pd.DataFrame()
i=1
for ch in tqdm(range(len(all_chrs))):
    #print(i)
    src = get_segs(atac1,all_chrs[ch])
    tmp = clean_segments(src)
    for t in tmp:
        results = results.append({'chr' : all_chrs[ch], 'start' : t[0], 'stop' : t[1]}, ignore_index = True)
    i=i+1

results.to_csv("sp4.lncRNAs_overlap_processed.csv")
