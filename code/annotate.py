import pandas as pd
from collections import Counter

#PROCESSING DATA

#process_raw and data files need to be in the folder you are running this script in

import process_raw #import process_raw.py, also on github

loci = pd.read_csv('top_fst_2pops_loci')
df = pd.read_csv('all_annotations.txt', sep="\t", header=None, names=["chr","source","type","start","stop","s1","s2","s3","info"])
df.drop(df[df['type'] == "region"].index, inplace = True)
df['info'] = df['info'].str.replace(',',';') #will help with parsing region info

# first column of loci df needs to be chr_pos for program to work
#loci["posi"]=loci["posi"].astype(int)
#loci["strpos"]=loci["posi"].astype(str)
#loci["ID"]=loci["chr"]+"_"+loci["strpos"]
#cols = list(loci.columns)
#cols = [cols[-1]] + cols[:-1]
#loci = loci[cols]
#loci["posi"]=loci["posi"].astype(str)

#ANNOTATING

annot_df=process_raw.annotate_raw(df,loci, one_length_print=False) #function defined in process_raw.py
annot_df.pos = annot_df.pos.astype(int)
annot_df.to_csv("annotated01.csv") #basic annotation, includes overlaps in one row
annot_nooverl=annot_df[annot_df['region']!="genes_overlap"]
print("done basic annotation")

overl = annot_df[annot_df['region']=="genes_overlap"]
overlapping_loci=process_raw.process_overlap(df,overl)
cdf=pd.concat([annot_nooverl,overlapping_loci]) #for now we keep both genes if SNPs falls in both
cdf.pos = cdf.pos.astype(int)
cdf.to_csv("annotated02.csv") #basic annotation, overlaps are resolved into multiple rows
print("done processing overlaps annotation")

promoters=process_raw.process_annotation_data(df)
missing_annot = annot_df[annot_df['region']=="not_annot"]
promoter_loci=process_raw.process_promoter(missing_annot,promoters)
promoter_loci.pos = promoter_loci.pos.astype(int) # loci in promoters, includes promoter-promoter overlaps in one row
print("done promoter processing")

poverl = promoter_loci[promoter_loci['region']=="genes_overlap"]
poverlapping_loci=process_raw.process_overlap(promoters,poverl)
cdf2=pd.concat([promoter_loci,poverlapping_loci])
cdf2.pos = cdf2.pos.astype(int)
cdf2.to_csv("promoters02.csv") # loci in promoters, promoter-promoter overlaps are resolved into mutliple rows
print("done processing overlaps in promoter annotation")

# we accounted for pos that fall in more than 1 gene region, and pos that fall in more than 1 promoter region
# now we correct for if pos falls into both gene body and promoter region

new_df=process_raw.annotate_raw(promoters,loci, one_length_print=True) 
new_df = new_df[new_df['region']=="genes_overlap"] #now it can be gene-gene, gene-promoter, gene-gene-promoter, gene-promoter-promoter or promoter-promoter overlap
new_df=new_df[~new_df.pos.isin(poverl.pos)] #remove promoter-promoter overlap
new_loci=process_raw.process_overlap(promoters,new_df)

ori_overlaps=annot_df[annot_df["region"]=="genes_overlap"]
oridf=cdf[cdf.pos.isin(ori_overlaps.pos)] # get annotations for overlapped genes

new_loci.gene = new_loci.gene.fillna('NaN')
oridf.gene = oridf.gene.fillna('NaN')
annot_df.gene = annot_df.gene.fillna('NaN')

new_loci.to_csv("new_loci.csv")
print(new_loci.head())

print("starting gene-promoter overlap resolve")
outdf=pd.DataFrame()
for index, row in new_loci.iterrows():
  found_in_oridf = oridf[(oridf["pos"]==row[2]) & (oridf["chr"]==row[0])] # chr, pos also in origianl df containing overlaps (no promoter)
  others_in_new_loci = new_loci[(new_loci["pos"]==row[2]) & (new_loci["chr"]==row[0])] # other same chr, pos in this df
  if len(found_in_oridf) > 0: # gene-gene overlap?
    if len(found_in_oridf) == len(others_in_new_loci): # regular gene-gene overlap, already in annotated02.csv
      pass
    elif len(found_in_oridf) < len(others_in_new_loci): #special case, gene-gene-promoter, resolve later if it comes up
      print("Error: more hits in new_loci than in oridf, due to gene-gene-promoter overlap")
    else:
      print("Error: more hits in oridf than in new_loci")
  else: # gene-promoter or gene-promoter-promoter, this is what we are really interested in
    dfmatch=annot_df[(annot_df["chr"]==row[0]) & (annot_df["pos"]==row[2]) & (annot_df["gene"]==row[1])] # chr, pos, gene name in original df
    if len(dfmatch) == 1: # this is the original gene name, hit for the gene body
      outdf = outdf.append({'chr' : row[0], 'pos' : row[2], 'region' : list(dfmatch.region)[0], "gene" : list(dfmatch.gene)[0]}, ignore_index = True)
    if len(dfmatch) == 0: # this is the new gene hit, hit for the promoter
        r="promoter" + ":" + row[3] # region it would fall in in "extended gene body", if intron, UTR, etc. if is promoter for regular gene
        outdf = outdf.append({'chr' : row[0], 'pos' : row[2], 'region' : r, "gene" : row[1]}, ignore_index = True)
    if len(dfmatch) > 1:
      print("Error: more than one hit in original df")

outdf[['promoter','p_type']] = outdf['region'].str.split(':',expand=True)

print("checking gene-promoter output df")
tdf=outdf.groupby(['chr','pos']).count()
tdf["ok"]= tdf["region"]-tdf["p_type"] # if at least 1 -> there is a nonpromoter region
print(len(tdf[tdf["ok"]<1])) # should be 0, if not, nonpromoter is missing somewhere
print(len(tdf[tdf["p_type"]<1])) # should be 0, if not, promoter is missing somewhere

p_overl_g=outdf[outdf["promoter"]=="promoter"]
p_overl_g["p_type"] = p_overl_g["p_type"].map({"intron":"protein_coding", "3'UTR":"protein_coding",
                                               "5'UTR":"protein_coding", "alternative UTR":"protein_coding",
                                               "exon":"protein_coding", "tRNA":"tRNA", "lnc_RNA":"lnc_RNA", "miRNA":"miRNA", "pseudogene":"pseudogene",
                                               "snRNA":"snRNA", "snoRNA":"snoRNA", "rRNA":"rRNA"})
p_overl_g.drop('region', axis=1, inplace=True)
p_overl_g.drop('promoter', axis=1, inplace=True)
p_overl_g["region"]=p_overl_g["p_type"]
p_overl_g.drop('p_type', axis=1, inplace=True)
p_overl_g.to_csv("promoter_gene_overl.csv") # new promoter regions I didn't find before because they are gene-promoter-promoter or gene-promoter overlaps

#CHECK POINT
print("starting check points")
only_promoter_loci=promoter_loci[promoter_loci['region']!="not_annot"] #loci in promoters
only_not_annot_prom=promoter_loci[promoter_loci['region']=="not_annot"] #loci not annotated even after promoter processing
only_annot_new_df=annot_df[annot_df['region']!="not_annot"] #loci anywhere else

check1=len(only_promoter_loci) + len(only_annot_new_df) + len(only_not_annot_prom)
check2=len(loci)
print("is everything in order?")
print(check1)
print(check2)

only_promoter_loci.to_csv("promoters01.csv") #includes all SNP that fall in promoters
print("done done")
