import pandas as pd
import numpy as np

add_promoters=True

df = pd.read_csv('all_annotations.txt', sep="\t", header=None, names=["chr","source","type","start","stop","s1","s2","s3","info"], skiprows=0)

print(df.head())

df.drop(df[df['type'] == "region"].index, inplace = True)
imp_list=["gene","snRNA", "tRNA", "lnc_RNA", "snoRNA", "rRNA", "miRNA"]

df=df[df["type"].isin(imp_list)] # these are the ones for which we looked at promoters for

if add_promoters:
    #adding promoters
    for l in imp_list:
        df.loc[(df["type"]==l) & (df["s2"]=="+"), "start"] -= 2000
        df.loc[(df["type"]==l) & (df["s2"]=="-"), "stop"] += 2000

    #fix promoters
    # as sometimes when I add and substract above, I get off the chromosome, either get a negative number or run over the end of the chromosome
    # first, fix going over max
    chr_lens=pd.read_csv('chr_lengths', sep=" ", header=None, names=["n","chr","start","stop"])

    chrlens_dic={}
    for index, row in chr_lens.iterrows():
        chrlens_dic[row[1]]=row[3]

    df["temp"]=df["chr"]
    df["temp"]=df['temp'].map(chrlens_dic)
    df["maxstop"]=df["temp"] - df["stop"]
    df["newstop"]=df["stop"]
    df.newstop = np.where(df.maxstop < 0, df.temp, df.stop)
    df["stop"]=df["newstop"]
    df=df.drop(columns=['temp','maxstop','newstop'])

    #second, fix going under
    df["newstart"]=df["start"]
    df.newstart = np.where(df.newstart < 0, 0, df.newstart)
    df["start"]=df["newstart"]
    df=df.drop(columns=['newstart'])


df=df.drop(columns=['source','type','s1', 's2', 's3','info'])

df=df.sort_values(by=["chr",'start'])
all_chrs=df.chr.unique()

#check if there is overlap between regions of the same set

def clean_segments(src):
    if len(src) <= 1:
        return src

    a = src[0]
    dest = []
    for b in src[1:]:
        # print(a,b, check_overlap(a,b))
        if check_overlap(a,b):
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
for ch in all_chrs:
    #print(i)
    src = get_segs(df,ch)
    tmp = clean_segments(src)
    for t in tmp:
        results = results.append({'chr' : ch, 'start' : t[0], 'stop' : t[1]}, ignore_index = True)
    i=i+1


results["length"]=results["stop"]-results["start"]
print(results.length.sum())
