# From raw reads to SNPs

## Checking quality of sequencing data

```
pip install multiqc
spack load fastqc@0.11.7

-----------
#!/bin/bash

yourfilenames=`ls /users/c/p/cpetak/WGS/all_fastqs/18170X*.fastq`

for file in $yourfilenames

do
	fastqc $file -o /users/c/p/cpetak/WGS/fastqc_output/
done
-----------

cd /users/c/p/cpetak/WGS/fastqc_output
multiqc .
```

### Results:

[Multiqc Report nicely displayed is available here](https://htmlpreview.github.io/?https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/data/multiqc_report.html) 

## Mapping to the reference genome

```
spack load bwa@0.7.17
spack load samtools@1.10
bwa index GCF_000002235.5_Spur_5.0_genomic.fna

-----------
while read line ; do
        F1=$(cut -d ' ' -f1 <<< $line)
        F2=$(cut -d ' ' -f2 <<< $line)
        echo "$F1 -- $F2"
        FILE=$(mktemp)
        cat header.txt >> $FILE
        echo "spack load samtools@1.10" >> $FILE
        echo "spack load bwa@0.7.17" >> $FILE
        ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna"
        out_name=$(cut -d '.' -f1 <<< $F1)
        echo "bwa mem -t 1 -M $ref /users/c/p/cpetak/WGS/all_fastqs/$F1 /users/c/p/cpetak/WGS/all_fastqs/$F2 | samtools view -S -b > /users/c/p/cpetak/WGS/BWA_out/$out_name.bam" >> $FILE
          sbatch $FILE
          sleep 0.5
          rm $FILE
done < $1
-----------
```

## Checking mapping statistics

```
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools sort /users/c/p/cpetak/WGS/BWA_out/$line -o /users/c/p/cpetak/WGS/BWA_out/$out_name.sorted.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools rmdup /users/c/p/cpetak/WGS/BWA_out/$line /users/c/p/cpetak/WGS/BWA_out/$out_name.rmdup.bam" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools flagstat /users/c/p/cpetak/WGS/BWA_out/$line | awk 'NR>=6&&NR<=13 {print \$1}' | column -x >> /users/c/p/cpetak/WGS/$out_name.flagstats.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
-----------
while read line ; do
	echo "$line"
	FILE=$(mktemp)
  	cat header.txt >> $FILE
	echo "spack load samtools@1.10" >> $FILE
	out_name=$(cut -d '.' -f1 <<< $line)
	echo "samtools depth /users/c/p/cpetak/WGS/BWA_out/$line | awk '{sum+=\$3} END {print sum/NR}' >> /users/c/p/cpetak/WGS/$out_name.coverage.txt" >> $FILE
  	sbatch $FILE
  	sleep 0.5
  	rm $FILE
done < $1
-----------
```

### Results:

In all 3 images below, x axis is the 140 individuals

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/coverage_fig.png" width="400" />

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/flagstat_fig.png" width="400" />

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/mapping_stat.png" width="400" />

## ANGSD for PCA

```bash
cd /users/c/p/cpetak/WGS/angsd

ref="/users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna" 
# latest version of the reference genome downloaded from NCBI on the 7th of October 2020
# angsd version: 0.933-102-g7d57642 (htslib: 1.11-9-g2264113) build(Oct 16 2020 18:14:45)

./angsd -b /users/c/p/cpetak/WGS/all_rmdups_jo.txt \
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \ # 85% of all individuals (140)
-setMinDepthInd 4 \ # note that later we'll use 3 here, filtering is stricter for now to reduce data for PCA
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doGlf 2 \ # gives us the Beagle format which will be used by pcangsd
-SNP_pval 1e-6
	
```

```bash
python /users/c/p/cpetak/pcangsd/pcangsd.py -beagle /users/c/p/cpetak/WGS/allpopstrict_angsd_polysites.beagle.gz -o /users/c/p/cpetak/WGS/pcangsd_covmatrix -threads 16
```

```R
#R
C <- as.matrix(read.table("pcangsd_covmatrix.cov"))
ids <- read.table("~/Downloads/pca_pops.txt") #text file with 20 lines of the single word BOD, then 20 lines of CAP etc in the order they appeared in all_rmdups_jo.txt
e <- eigen(C)
# base R
plot(e$vectors[,1:2],xlab="PC1",ylab="PC2", bg=ids$V1, pch=21)
#ggplot
library(ggplot2)
library(tidyverse)
df <- data.frame(pop = ids$V1, PC1 = e$vectors[,1], PC2 = e$vectors[,2])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()
```

### Results:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/PCA_1.png" width="400" />
	
3 individuals seem to be very different from the other 137 individuals. Thus, these 3 were dropped from further analysis. New PCA with 137 individuals (ANGSD was rerun with only 137 individuals):

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/PCA_2.png" width="400" />

No clustering by population can be seen.

### PCA for quality check

Histogram of average coverage per individual:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/hist_coverage.png" width="400" />
	

```
#R
covdata <- read.table("covs.txt") #this is without outliers, average coverage per individual
covdata$V2 <- ifelse(covdata$V1<6.5, "little", "lot")

ids <-covdata

#ggplot
df <- data.frame(pop = ids$V2, PC1 = e$vectors[,1], PC2 = e$vectors[,2])
df= rownames_to_column(df)
ggplot(df, aes(x=PC1, y=PC2, fill=pop)) +
  geom_point(size=3, shape=21) +
  theme_bw()
```

PCA of average coverage:

<img src="https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/images/PCA_cov.png" width="400" />
     
Again, no clustering is visible.

## ANGSD for SNPs

First, run ANGSD to get bcf file with per site information for each individual from all populations.


```bash
./angsd -b /users/c/p/cpetak/WGS/make_vcf/list_for_angsd.txt \ 
-ref ${ref} \
-anc ${ref} \
-out /users/c/p/cpetak/WGS/make_vcf/all_pop_angsd \
-nThreads 16 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 30 \
-minQ 20 \
-minInd 119 \
-setMinDepthInd 3 \
-skipTriallelic 1 \
-dobcf 1 \
-GL 1 \
-doPost 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
-SNP_pval 1e-6
```

-> output of interest: all_pop_angsd.bcf

```bash
# convert to vcf
spack load bcftools@1.10.2
bcftools view all_pop_angsd.bcf > all_pop_angsd.vcf
# then to add GT info
vcfglxgt all_pop_angsd.vcf > fixed_all_pop_angsd.vcf
# then to account for missing information
sed -i 's/0\/0\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/NA\:0\,0\,0\:0\,0\,0\:0\,0\,0\:0/g' fixed_all_pop_angsd_copy.vcf
# then keep only GT information
awk -v FS="\t" -v OFS="\t" '{for(i=9;i<=NF;i++) {split($i, gt, ":"); $i=gt[1]} print}' fixed_all_pop_angsd_copy.vcf > fixed_all_pop_angsd_copy_onlyGT.vcf
# getting rid of vcf header
head -916 fixed_all_pop_angsd_copy_onlyGT.vcf > vcf_head.vcf
cat fixed_all_pop_angsd_copy_onlyGT.vcf | grep -v "#" > vcf_tail.vcf 
```

### Filtering for biallelic

```bash
# getting first column of vcf as chromosome.position
awk -F "\t" '{print $1$2, $0}' vcf_tail.vcf > vcf_tail_idd2
# keeping only lines in vcf that are also in list of filtered positions
awk 'FNR==NR{a[$0];next}($1 in a)' filt_posi_fixed* vcf_tail_idd2 > filtered_vcf
# where filt_posi_fixed is the list of loci after filtering, 1 column, chromosome.position format
# results in a file that has 991,430 positions (994,220 only MAF filter), because remember here we also filtered for SNP_pval
```

*This is coming from the filtering steps described in [Code for Filtering steps](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/Filtering_steps.md)

## Site frequency spectrums and Fst

First, I run angsd on each population separately. E.g.

```bash
./angsd -b /users/c/p/cpetak/WGS/BOD_rmdups_jo.txt 
-ref /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-anc /users/c/p/cpetak/WGS/reference_genome/GCF_000002235.5_Spur_5.0_genomic.fna 
-out /users/c/p/cpetak/WGS/angsd_new/BOD_angsd_allsites 
-nThreads 16 
-remove_bads 1 
-C 50 
-baq 1 
-minMapQ 30 
-minQ 20 
-minInd 17 # 85% of 20 individuals
-setMinDepthInd 3 
-skipTriallelic 1 
-GL 1 
-doCounts 1 
-doMajorMinor 1 
-doMaf 1 
-doSaf 1 
-doHWE 1
```

For FOG and CAP, I also run the above code without outliers 

Since I have 7 populations, I have 21 possible pairs of populations. For each of the possible pairs:

```bash
#!/bin/sh
dir=/users/c/p/cpetak/WGS/angsd_new
while read line ; do #give this script a list of pop pairs displayed as pop1.pop2
    pop1=$(cut -d '.' -f1 <<< $line)
    pop2=$(cut -d '.' -f2 <<< $line)
    echo $pop1
    echo $pop2
    FILE=$(mktemp)
    cat header.txt >> $FILE
    echo "cd /users/c/p/cpetak/WGS/angsd/misc" >> $FILE
    echo "./realSFS ${dir}/${pop1}_angsd_allsites.saf.idx ${dir}/${pop2}_angsd_allsites.saf.idx -P 16 -fold 1 > ${dir}/pairwise_fst/${pop1}_${pop2}_allsites.sfs" >> $FILE #folded option!
    sbatch $FILE
    sleep 0.5
    rm $FILE
done < $1
```

For all pairs containing CAP and FOG the angsd output with no outliers was used (from angsd_noout)

Then for each pair (using TER.BOD as an example)

```bash
cd /users/c/p/cpetak/WGS/angsd/misc
dir=/users/c/p/cpetak/WGS/angsd_new
./realSFS fst index ${dir}/TER_angsd_allsites.saf.idx ${dir}/BOD_angsd_allsites.saf.idx -sfs ${dir}/pairwise_fst/TER_BOD_allsites.sfs -fold 1 -fstout ${dir}/pairwise_fst/TER_BOD_allsites -whichFst 1
```

Finally

```bash
./realSFS fst print /users/c/p/cpetak/WGS/angsd_new/pairwise_fst/TER_BOD_allsites.fst.idx > /users/c/p/cpetak/WGS/angsd_new/pairwise_fst/TER_BOD_allsites.fst
```

To filter by MAF 

```bash
sed s/"\.1"/"\.1,"/g bayenv_onlypos_0025filter.csv* > bayenv_onlypos_temp.csv
awk -F "," '{print $2"_"$1, $0}' bayenv_onlypos_temp.csv > bayenv_onlypos_temp2.csv
cut -d' ' -f1 bayenv_onlypos_temp2.csv > bayenv_onlypos_temp3.csv
# now it is one column, with pos_chr format, rm intermediate temp files and rename to _cleaned
# repeat with each .fst file too (not exactly as they have different starting format) to match format:
awk '{print $2"_"$1, $0}' CAP_BOD_allsites.fst # example file
# then sort each of the fst files along with the bayenv file
sort -k1 
# for each fst file
join -1 1 -2 1 sorted_bayenv_onlypos_cleaned.csv ${pop1}_${pop2}_allsites_cleaned.fst > joined_${pop1}_${pop2}_allsites
# IMPORTANT! Make sure that the column you are joining on (in this case the first column) has the same format in both files you are joining! E.g. pos_chr. I chose this IDing instead of chr_pos to avoid sorting issues.
```

*This is coming from the filtering steps described in [Code for Filtering steps](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/Filtering_steps.md)

### Global pairwise Fst

```bash
./realSFS fst stats ${dir2}/pairwise_fst_cleaned/${pop1}_${pop2}_allsites.fst.idx > ${dir2}/pairwise_fst_cleaned/${pop1}_${pop2}_global.fst
```

-> output: unweighted and weighted (second preferred) global Fst for each pair. 

### Per-site nucleotide diversity

Calculating sfs for each pop, reusing code from doing it per population pair.

```bash
./realSFS ${dir}/${line}_angsd_allsites.saf.idx -P 16 -fold 1 > ${dir}/${line}_allsites.sfs
```

Then, calling saf2theta

```bash
./realSFS saf2theta ${dir}/${line}_angsd_allsites.saf.idx -sfs ${dir}/${line}_allsites.sfs -outname ${dir}/${line}_thetaout_div
```

Then, to calculate Tajimas D

```bash
./thetaStat do_stat ${dir}/${line}_thetaout_div.thetas.idx
.thetaStat do_stat out.thetas.idx -win 100 -step 100  -outnames theta.thetasWindow.gz
```

then repeated but with -win 1000 -step 1000 -> eg BOD_1000.thetasWindow.gz.pestPG

### Linkage disequilibrium

```bash
# PLINK 1.9
# PLINK v1.90b6.25 64-bit (5 Mar 2022)
# cp filtered_vcf, rename as temp_vcf
cut -d" " -f2 temp_vcf > proc_temp_vcf
sed -i 's/_//g' vcf_header_temp
tail -10000 proc_temp_vcf > temp_vcf2 # downsample as it is too big
cat vcf_header_temp temp_vcf2 > temp_vcf_ready
sed -i 's/NA/./g' temp_vcf_ready

# LD decay
head -500 temp_vcf_ready > vcf_input # they were all on the same chromosome, NW_022145612.1
./plink --vcf vcf_input --recode --allow-extra-chr --out temp_vcf_plink #now we have pad and map
./plink --file temp_vcf_plink --make-bed --allow-extra-chr --out afterQC
# adjust the number of SNPs and inter-SNP distances for which you want to compute LD
./plink --bfile afterQC --r2 --allow-extra-chr --ld-window-r2 0 --ld-window 10000 --ld-window-kb 500000 --out resultLD

# Create heatmap
./plink --bfile afterQC --r2 square --allow-extra-chr --out resultLD_heat
awk -F "\t" '{print $2}' vcf_input > locs_info

# LD pruning
cat vcf_header_temp proc_temp_vcf > temp_vcf_ready # no downsampling this time
sed -i 's/NA/./g' temp_vcf_ready
./plink --vcf temp_vcf_ready --recode --allow-extra-chr --out temp_vcf_plink
./plink --file temp_vcf_plink --make-bed --allow-extra-chr --set-missing-var-ids @:#[b37]\$1,\$2 --out afterQC # my SNPs didn't have IDs
./plink --bfile afterQC --allow-extra-chr --indep-pairwise 50 5 0.5 
```

