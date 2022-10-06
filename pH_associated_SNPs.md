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

OutFlank fit:

2 pops:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/outflank_2pop_default.jpg" width="400" />

Distribution of p-values:

2 pops:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/p_hist_2pop_default.jpg" width="400" />

note the excess of p-values near 1, indicating poor fit of the left tail. Unfortunately, this is a fault of OutFlank, "OutFLANK will not fit the left tail of the FST distribution well".

For the 2 populations one: 2241 outliers with default settings  (q, Hmin, trims, everything) 

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/outliers_2pops_default.jpg" width="400" />

## Analysis of outlier loci

separate?
