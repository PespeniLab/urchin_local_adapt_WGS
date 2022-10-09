# Gether all annotation data (from version 5.0) into 1 dataframe

import pandas as pd
import numpy as np

# Make master df
df=pd.read_csv("annotated02.csv")
df = df.iloc[: , 1:]
df = df[["chr", "pos", "region","gene"]]
df=df.groupby(["chr","pos"]).agg(lambda x: list(x))

p2=pd.read_csv("promoters02.csv")
p2 = p2.iloc[: , 1:]

# Not annotated
not_annot=p2[p2["region"]=="not_annot"]
df=df.merge(not_annot, how='outer', on=['chr', 'pos'])
df.drop('gene_y', axis=1, inplace=True)
df=df.rename(columns={"region_y": "Not_annot"})
df.Not_annot = df.Not_annot.fillna('NaN')
df["Not_annot"] = df["Not_annot"].map({"not_annot":1, "NaN":0,})

# Promoters
proms=p2[p2["region"]!="not_annot"]
proms=proms[proms["region"]!="genes_overlap"]
proms["region"] = proms["region"].map({"protein_coding":"protein_coding", "lnc_RNA":"lnc_RNA","tRNA":"tRNA","intron":"protein_coding" })
proms2=proms.groupby(["chr","pos"]).agg(lambda x: list(x))

df=df.merge(proms2, how='outer', on=['chr', 'pos'])
df=df.rename(columns={"gene": "prom_gene", "region":"prom_region"})
#df['promoter'] = df["promoter"].where(df["promoter"].isnull(), 1).fillna(0).astype(int)

# Adding overlap info (gene-promoter, gene-promoter-promoter)
pgo=pd.read_csv("promoter_gene_overl.csv")
pgo = pgo.iloc[: , 1:]
pgo=pgo.groupby(["chr","pos"]).agg(lambda x: list(x))
df=df.merge(pgo, how='outer', on=['chr', 'pos'])
df=df.rename(columns={"gene": "prom_gene_overl", "region":"prom_region_overl"})
df['prom_gene']=df['prom_gene'].fillna(df['prom_gene_overl'])
df['prom_region']=df['prom_region'].fillna(df['prom_region_overl'])
df.drop('prom_gene_overl', axis=1, inplace=True)
df.drop('prom_region_overl', axis=1, inplace=True)

df["promoter"]=df["prom_region"]
df['promoter'] = df['promoter'].str.len()
df['promoter'] = df['promoter'].fillna(0)

possibilities=["pseudogene", "snRNA", "tRNA", "lnc_RNA", "snoRNA", "rRNA", "miRNA", "exon", "3'UTR", "5'UTR", "alternative UTR", "intron"]
for p in possibilities:
  df[p]=df.region_x.map(set([p]).issubset)
  df[p] = df[p].astype(int)

# account for cases where gene-gene overlap happened between the same regions

df["num_gene_hit"]=df.region_x.str.len()
df["num_gene_hit"] = df["num_gene_hit"] - df["Not_annot"]
df["accunted_for"] = df["pseudogene"] + df["snRNA"] + df["tRNA"] + df["lnc_RNA"] + df["snoRNA"] + df["rRNA"] + df["miRNA"] + df["exon"] + df["3'UTR"] + df["5'UTR"] + df["alternative UTR"] + df["intron"]
df["check"]=df["num_gene_hit"] - df["accunted_for"]
cdf=df[(df["check"]!=0) & (df["promoter"]==0)]

from collections import Counter
for index, row in cdf.iterrows():
    c=Counter(row[2])
    for e in row[2]:
        df.loc[index,e]=c[e]

print(len(df[(df["num_gene_hit"]==0) & (df["accunted_for"]==0) & (df["Not_annot"]==0)])) #should be 0
df["accunted_for"] = df["pseudogene"] + df["snRNA"] + df["tRNA"] + df["lnc_RNA"] + df["snoRNA"] + df["rRNA"] + df["miRNA"] + df["exon"] + df["3'UTR"] + df["5'UTR"] + df["alternative UTR"] + df["intron"]
df["check"]=df["num_gene_hit"] - df["accunted_for"]
print(len(df[(df["check"]!=0) & (df["promoter"]==0)])) #should be 0

df.to_csv("gathered_annotation.csv")
tdf=df[df["Not_annot"]==1]
tdf.to_csv("input_for_extra_regannot_notannot.csv") # to be translated into the 3.1 version

# now if I do sum across a column, I get the total number of eg intron hits (so more than one if a pos fell in intron-inton overlap)
# if instead I do intron != 0 I get number of pos that fell into an intron
# similarly, if I do sum across promoters I get total number of promoter hits
# but if I do promoter != 0 I get number of pos that fell into promoter (tho they could have also fell into an intron for example)
