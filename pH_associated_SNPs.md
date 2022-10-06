# Identification of pH variability associated SNPs

## pH data analysis

See notebook code/ph_data_clean.ipynb

## Outflank

Cleaning filtered_vcf (output of ANGSD) file from Part 1:

```bash
#cleaning for outflank
sed -i 's/NA/NA\/NA/g' filtered_vcf
cut --complement -d$' ' -f1 filtered_vcf > cleaned_filtered_vcf
#taking out outliers
cut --complement -d$'\t' -f30,56,57 cleaned_filtered_vcf > cut_filtered_vcf
# put back vcf header
cat vcf_head.vcf cut_filtered_vcf > topped_vcf.vcf
# take out names of outliers in vcf header
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/CAP\_18170X101\_200925\_A00421\_0244\_AHKML5DSXY\_S121\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/FOG\_18170X127\_200925\_A00421\_0244\_AHKML5DSXY\_S147\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
sed -i 's/\/users\/c\/p\/cpetak\/WGS\/BWA\_out\/FOG\_18170X128\_200925\_A00421\_0244\_AHKML5DSXY\_S148\_L002\_R1\_001\.rmdup\.bam//g' topped_vcf.vcf
# removing accidental double tabs
sed 's:\t\t*:\t:g' topped_vcf.vcf > test.vcf
```

Outflank

```R
#calculating Fst values
library(OutFLANK)
library(vcfR)

vcf <- read.vcfR("test.vcf", verbose=FALSE)
ind <- read.table("Pop.txt", header=TRUE)

convertVCFtoCount3 <- function(string){
    a <- as.numeric(unlist(strsplit(string, split = c("[|///]"))))
    odd = seq(1, length(a), by=2)
    a[odd] + a[odd+1]
}
all.vcf.gen <- vcf@gt[,-1]
system.time(gen_table <- matrix(convertVCFtoCount3(all.vcf.gen), ncol=ncol(all.vcf.gen)))

locinames <- paste(vcf@fix[,"CHROM"], vcf@fix[,"POS"], sep="_")
SNPdata <- t(gen_table)
SNPdata[is.na(SNPdata)] <- 9
k <- max(ind$pop)

FstDataFrame <- MakeDiploidFSTMat(SNPdata,locinames,ind$pop)

write.csv(FstDataFrame, file = "results.csv")
```

Where Pop.txt is just a list of 57 1s (first "pop" BOD + CAP + FOG individuals) and 80 2s (second "pop").

```R
#quality checking
library(OutFLANK)
library(vcfR)

FstDataFrame <- read.csv(file = 'results.csv', header=TRUE,row.names=1)

pdf("line.pdf")
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, xlim=c(-0.01,0.3), ylim=c(-0.01,0.3), pch=20)
abline(0,1)
dev.off()

pdf("dots.pdf")
plot(FstDataFrame$He, FstDataFrame$FSTNoCorr, pch=20, col="grey")
dev.off()

pdf("hist.pdf")
hist(FstDataFrame$FSTNoCorr[FstDataFrame$He>0.1],xlim=c(0,0.3), breaks=50)
dev.off()
```

```R
#searching for outliers
FstDataFrame <- read.csv(file = 'results.csv', header=TRUE,row.names=1)
k=2
q=0.05
outlier <- OutFLANK(FstDataFrame, NumberOfSamples=k) #investigate options
write.csv(outlier, file = "outliers_from_outflank.csv")

pdf("outflank.pdf")
OutFLANKResultsPlotter(outlier, withOutliers = TRUE,NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =FALSE, RightZoomFraction = 0.05, titletext = NULL) #investigate options
dev.off()

pdf("p_hist.pdf")
hist(outlier$results$pvaluesRightTail)
dev.off()

num_out <- sum(outlier$results$qvalues<q, na.rm=TRUE)

if (num_out > 0) {
print("there are outliers:")
print(num_out)
pdf("outliers.pdf")
plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
    points(outlier$results$He[outlier$results$qvalues<q], outlier$results$FST[outlier$results$qvalues<q], pch=21, col="blue")
dev.off()

top_candidates <- outlier$results$qvalues<q & outlier$results$He>0.1
topcan <- outlier$results[top_candidates,]

write.csv(topcan, file = "top_fst.csv")
}
```

Line and dot plots looked good. A description of what these plots show: http://rstudio-pubs-static.s3.amazonaws.com/305384_9aee1c1046394fb9bd8e449453d72847.html

Line plot

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/line_2pop.jpg" width="400" />

Dot plot

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/dots_2pop.jpg" width="400" />

This is just meant to show how for low He SNPs Fst is inflated. SNPs for which He < 0.1 were removed in outlier calculations (OutFlank function, default). Note how there aren't any positions He < 0.05, as there sites were filtered out already.

Now let's look at the distribution of Fsts (He > 0.1):

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/hist_2pop.jpg" width="400" />

Running OutFlank of default settings (as written above in code chunk) resulted in a relatively good fit for both of the distributions.

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/outflank_2pop_default.jpg" width="400" />

Distribution of p-values:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/p_hist_2pop_default.jpg" width="400" />

note the excess of p-values near 1, indicating poor fit of the left tail. Unfortunately, this is a fault of OutFlank, "OutFLANK will not fit the left tail of the FST distribution well".

For the 2 populations one: 2241 outliers with default settings  (q, Hmin, trims, everything) 

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/outliers_2pops_default.jpg" width="400" />



## Pairwise Fst

After calculating an Fst value using realsfs (see Part 1 of analysis) for each locus for each pair of population, I created 2 files - one with the list of pop pairs within pH groups (i.e. low pH variability pop vs low pH variability pop or high pH variability pop vs high pH variability pop), and one with a list of pop pairs between pH groups (i.e. low pH variability pop vs high pH variability pop). For each of these categories I combined the corresponding joined files. 

```bash
cat $(cat between_pairs) > combined_between
cat $(cat within_pairs) > combined_within
```

 I then averaged Fsts coming from within and between pop pairs separately to following way:

```python
import pandas as pd
import matplotlib.pylab as plt
import numpy as np

col_names1=["ID","chr","pos","A", "B"]

df = pd.read_csv('combined_between',names=col_names1, sep="\s")
df["pos"]=df.pos.astype('int64', copy=False)

df["fst"]=df["A"]/df["B"] # see realsfs documentation for the meaning of these values
df.replace([np.inf, -np.inf], np.nan, inplace=True)
df['fst'] = df['fst'].fillna(0)

df = df.drop(["ID","A","B"], 1)

ndf=df.groupby(["chr","pos"]).mean()

ndf.to_csv('average_between.csv') # do same with within
```

Sometimes the average Fst is slightly negative, let's round that to 0

```bash
awk -F, '$3+0<0{$3=0}1' average_between.csv | sed 's/ /,/g' > average_between_cropped.csv
```

Next, we need to subtract average_within from average_between. I used the following code:

```python
import pandas as pd
import numpy as np

col_names1=["chr","pos", "fst"]

betw = pd.read_csv('average_between_cropped.csv', sep=",")
withi = pd.read_csv('average_within_cropped.csv', sep=",")

betw["pos"]=betw.pos.astype('int64', copy=False)
withi["pos"]=withi.pos.astype('int64', copy=False)
betw["fst"] = pd.to_numeric(betw["fst"], downcast="float")
withi["fst"] = pd.to_numeric(withi["fst"], downcast="float")

cdf=betw.merge(withi, on=["chr","pos"], how="outer")
cdf = cdf.fillna(0)

cdf["final_fst"]=cdf["fst_x"]-cdf["fst_y"]

cdf.to_csv("subbed.csv")
```

subbed.csv contains the final fst value for each locus.

Distribution of final pairwise Fst values:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/pairwise_FST_histogram.jpg" width="400" />

Then bootstrapped the top 1% cutoff point using the following code:

```python
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv("subbed.csv")

fsts=df.final_fst.to_numpy()

def bootstrap_sample(amounts):
    return np.random.choice(amounts, len(amounts), replace=True)

def percentile_99(sample):
     return np.percentile(sample, 99)

def bootstrap_confidence_interval(data):
    """
    Creates list of 10000 99th percentile bootstrap replicates.
    """

    itnum=10000

    bs_samples = np.empty(itnum)

    for i in range(itnum):
        bs_samples[i] = percentile_99(bootstrap_sample(data))

    return bs_samples

transactions_ci = bootstrap_confidence_interval(fsts)

print(np.percentile(transactions_ci, 95))

np.savetxt('99s.out', transactions_ci, delimiter=',')

fig3 = plt.figure()
plt.hist(transactions_ci, bins=40)
fig3.savefig('boots.png')
```

Output: 0.0260773276575

```bash
awk -F "," ' $6 >= 0.0260773276575 ' subbed.csv > pair_fst_outs
```

Number of outliers here: 9,780







## Analysis of outlier loci

separate?
