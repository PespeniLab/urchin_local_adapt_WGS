from scipy.stats import chisquare
import pandas as pd
import numpy as np

def do_chi(all_hit, tot):
  #calculated based on https://www.ncbi.nlm.nih.gov/genome/annotation_euk/Strongylocentrotus_purpuratus/102/
  #count x mean length (bp)
  bp_exon=245575 * 335
  bp_3UTR=40771628
  bp_5UTR=10192907
  bp_intron=218005 * 1808
  bp_lnc_paper=17703544 # number of basepairs (no double counting due to overlaps) from paper
  bp_lnc_ncbi=4850968 #only from ncbi
  bp_lnc=bp_lnc_paper + bp_lnc_ncbi
  bp_enhancer=5185974
  bp_promoter=51178642
  bp_noncoding=921855793-bp_exon-bp_3UTR-bp_5UTR-bp_intron-bp_lnc-bp_enhancer-bp_promoter
  all_len=[bp_exon,bp_3UTR,bp_5UTR,bp_intron,bp_lnc,bp_enhancer,bp_promoter,bp_noncoding]

  all_percent_expected=[]
  for a in all_len:
    all_percent_expected.append(a/921855793) #percentage of genome corresponding to each region

  #number of expected loci in each region if the high Fst loci were randomly distributed between the regions
  all_expected_rounded=[round(element * tot) for element in all_percent_expected]
  all_expected=[element * tot for element in all_percent_expected]

  print(chisquare(all_hit, f_exp=all_expected))
  print("[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter,num_noncoding]")
  print(all_hit)
  print(all_expected_rounded)
  print([i / j for i, j in zip(all_hit, all_expected)])

df=pd.read_csv("gathered_annotation.csv")
additional_lnc=pd.read_csv("lncdf.csv")
notannot_enhancers=pd.read_csv("reg_regions_withconfidence.csv")

#checking, should be of len 1
print(len(additional_lnc.groupby(by=["source_id","source_start"]).count().region.unique()))
print(len(notannot_enhancers.groupby(by=["chr50","pos50"]).count().region.unique()))

## Method 1: overlaps count

num_exon=len(df[df["exon"]!=0]) # could have overlapped with something else as well...
num_3UTR=len(df[df["3'UTR"]!=0])
num_5UTR=len(df[df["5'UTR"]!=0])
num_intron=len(df[df["intron"]!=0])
num_lnc=len(df[df["lnc_RNA"]!=0]) + len(additional_lnc)
num_enhancer=len(notannot_enhancers[notannot_enhancers["confidence"]>0]) #TODO
num_promoter=len(df[df["promoter"]!=0]) # could overlap with other promoter as well as gene

pseudogene=len(df[df["pseudogene"]!=0])
snRNA=len(df[df["snRNA"]!=0])
tRNA=len(df[df["tRNA"]!=0])
snoRNA=len(df[df["snoRNA"]!=0])
rRNA=len(df[df["rRNA"]!=0])
miRNA=len(df[df["miRNA"]!=0])
alternative_UTR=len(df[df["alternative UTR"]!=0])

num_noncoding=len(df[df["Not_annot"]==1]) - len(notannot_enhancers) - len(additional_lnc) + pseudogene + snRNA + tRNA + snoRNA + rRNA + miRNA + alternative_UTR
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter,num_noncoding]

print(sum(all_hit))
print(len(df))
print((sum(all_hit)-len(df))/len(df)*100) # it is this percent of an increase due to my method of counting

do_chi(all_hit, len(df))

## Method 2: everything is counted once, more conservative
# gene-gene neither is counted, promoter-promoter counted as 1, gene-promoter and gene-promoter-promoter only gene is counted.

# Make alternative df
tdf = pd.read_csv("annotated01.csv") #with overlaps in 1 row, called genes_overlap
tdf = tdf.iloc[: , 1:]
tdf = tdf[["chr", "pos", "region","gene"]]
print(len(tdf[tdf["region"]=="genes_overlap"])) #number of gene-gene overlaps
tdf = tdf[tdf["region"]!="genes_overlap"] # get rid of gene overlaps
tdf = tdf[tdf["region"]!="not_annot"] # keep only annotated

p=pd.read_csv("promoters02.csv") #with overlaps in 1 row, called genes_overlap
p = p.iloc[: , 1:]
notannot=p[p["region"]=="not_annot"] #keep only not annot

proms=p[p["region"]!="not_annot"] #keep only promoters
pt=proms.drop_duplicates(subset=['chr','pos']) #get rid of promoter-promoter overlaps
pd.options.mode.chained_assignment = None
pt["region"]="promoter"
df2=pd.concat([tdf, pt, notannot ])

# Count nucleotides in each region

num_exon=len(df2[df2["region"]=="exon"])
num_3UTR=len(df2[df2["region"]=="3'UTR"])
num_5UTR=len(df2[df2["region"]=="5'UTR"])
num_intron=len(df2[df2["region"]=="intron"])
num_lnc=len(df2[df2["region"]=="lnc_RNA"]) + len(additional_lnc)
num_enhancer=len(notannot_enhancers[notannot_enhancers["confidence"]>0]) #TODO
num_promoter=len(df2[df2["region"]=="promoter"])
num_noncoding=len(df) - sum([num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter])
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter,num_noncoding]

print(sum(all_hit))
print(len(df))
print((sum(all_hit)-len(df))/len(df)*100) # it is this percent of an increase due to my method of counting

do_chi(all_hit,len(df))

## Method 3: same as method 2 but for gene-gene overlap select 1 at random instead of discarding all
# Technically it drops first duplicate but since they were added to the csv randomly it shouldn't matter, so practically we end up selecting randomly

# Make alternative df
tdf = pd.read_csv("annotated02.csv") #with overlaps in multiple rows
tdf = tdf.iloc[: , 1:]
tdf = tdf[["chr", "pos", "region","gene"]]
print(len(tdf[tdf["region"]=="genes_overlap"])) #number of gene-gene overlaps
tdf=tdf.drop_duplicates(subset=['chr','pos']) # keep only first occurrence
tdf = tdf[tdf["region"]!="not_annot"] # keep only annotated

p=pd.read_csv("promoters02.csv") #with overlaps in 1 row, called genes_overlap
p = p.iloc[: , 1:]
notannot=p[p["region"]=="not_annot"] #keep only not annot

proms=p[p["region"]!="not_annot"] #keep only promoters
pt=proms.drop_duplicates(subset=['chr','pos']) #get rid of promoter-promoter overlaps
pd.options.mode.chained_assignment = None
pt["region"]="promoter"
df2=pd.concat([tdf, pt, notannot ])

# Count nucleotides in each region

num_exon=len(df2[df2["region"]=="exon"])
num_3UTR=len(df2[df2["region"]=="3'UTR"])
num_5UTR=len(df2[df2["region"]=="5'UTR"])
num_intron=len(df2[df2["region"]=="intron"])
num_lnc=len(df2[df2["region"]=="lnc_RNA"]) + len(additional_lnc)
num_enhancer=len(notannot_enhancers[notannot_enhancers["confidence"]>0]) #TODO
num_promoter=len(df2[df2["region"]=="promoter"])
num_noncoding=len(df) - sum([num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter])
all_hit=[num_exon,num_3UTR,num_5UTR,num_intron,num_lnc,num_enhancer,num_promoter,num_noncoding]

print(sum(all_hit))
print(len(df))
print((sum(all_hit)-len(df))/len(df)*100) # it is this percent of an increase due to my method of counting

do_chi(all_hit, len(df))
