import pandas as pd

#making this notebook to take out regions from possible enhancers
#this is NOT for checking Fst outliers in genes or checking if eg Fst is in UTR

df = pd.read_csv('all_annotations.txt', sep="\t", header=None, names=["chr","source","type","start","stop","s1","s2","s3","info"], skiprows=0)

#these are either overlaping or not important due to low sample size
df.drop(df[df['type'] == "region"].index, inplace = True)
df.drop(df[df['type'] == "mRNA"].index, inplace = True)
df.drop(df[df['type'] == "exon"].index, inplace = True)
df.drop(df[df['type'] == "CDS"].index, inplace = True)
df.drop(df[df['type'] == "sequence_feature"].index, inplace = True)
df.drop(df[df['type'] == "tandem_repeat"].index, inplace = True)
df.drop(df[df['type'] == "promoter"].index, inplace = True)
#so we are left with eg gene, lncRNA, rRNA snRNA, etc... we want to take out a region that overlaps with any of these as then they are not cis-regulatory regions

#to account for promoters, since we are searching for non promoter/intron enhancers
#if + strand -> substract 5000 from start
#else add 5000 to stop

df.loc[(df["type"]=="gene") & (df["s2"]=="+"), "start"] -= 5000

df.loc[(df["type"]=="gene") & (df["s2"]=="-"), "stop"] += 5000


df=df.drop(columns=['source','type','s1', 's2', 's3','info'])

#df["type"] = "annot"
df=df.sort_values(by=["chr",'start'])
all_chrs=df.chr.unique()
#df["length"] = df["stop"] - df["start"]

#CHECK IF THERE ARE NEGATIVE NUMBERS
#min(df["length"]) <= 0

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


results.to_csv("nonoverlapping_nopromoter_gffannotations.csv",index=False)
