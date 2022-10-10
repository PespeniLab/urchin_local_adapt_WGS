# From candidate SNPs to genomic regions and genes

## Chi-squared of genomic regions

#### Annotating outliers with v5.0 for nonregulatory variation

I annotated the outlier SNPs using the following code: [annotate.py](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/annotate.py), which uses [process_raw.](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/process_raw.py) These use the Spur v5.0 annotation available on NCBI RefSeq.

There are 6 types of annotations in various files as the output of above python code:

1. Gene body hits -> annotated01.csv, 
2. Hits in gene body - gene body overlaps are resolved into separate lines -> annotated02.csv
3. Promoter hits -> promoters01.csv, just defined as within 5000 bp of TSS.
4. Hits in promoter - promoter overlaps resolved into separate lines ->  promoters02.csv
5. Hits in gene body - promoter overlaps  and hits in gene body - promoter - promoter overlaps -> promoter_gene_overl.csv
6. Hits in gene body - gene body - promoter overlaps (rare, this annotation still shows pos as simple gene body - gene body overlap)

I gathered all of the above annotation data into one big dataframe using the [gather_annot_info.py](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/gather_annot_info.py) code. This created a dataframe with all outlier loci as rows and columns containing region information such as intron, promoter, etc with a number representing how many of these regions each loci "hit". So for example if a locus was found in 2 overlapping promoter regions, the whole row has 0s except for the promoter column which says 2.  [Here](https://github.com/Cpetak/urchin_adaptation/blob/main/data/2_pop_fst_step3/gathered_annotation.csv) is this file. From this dataframe, I then extracted rows that were not annotated at all -> input_for_extra_regannot_notannot.csv. Loci in this file were converted to v3.1 as follows. 

#### Annotating outliers with v3.1 for regulatory regions

Used this tool to convert locations to v3.1: https://www.ncbi.nlm.nih.gov/genome/tools/remap

Opened input_for_extra_regannot_notannot.csv in Visual Studio Code, and from there I copied into space in link above. Had to remove quatation marks and ".1" from the end of the chromosome numbers. Chr:pos column separated.

Regulatory regions from the literature: 

* ATAC-seq, DNA-seq: Shashikant, T., Khor, J. M., & Ettensohn, C. A. (2018). Global  analysis of primary mesenchyme cell cis-regulatory modules by chromatin  accessibility profiling. *BMC genomics*, *19*(1), 1-18.
* Chip-seq: Khor, J. M., Guerrero-Santoro, J., & Ettensohn, C. A. (2019).  Genome-wide identification of binding sites and gene targets of Alx1, a  pivotal regulator of echinoderm skeletogenesis. *Development*, *146*(16), dev180653.
* L. var. similarity: Tan, G., Polychronopoulos, D., & Lenhard, B. (2019). CNEr: A toolkit for exploring extreme noncoding conservation. *PLoS computational biology*, *15*(8), e1006940.
* Additional known enhancers: Arenas-Mena, C., Miljovska, S., Rice, E. J., Gurges, J., Shashikant, T., Wang, Z., ... & Danko, C. G. (2021). Identification and prediction  of developmental enhancers in sea urchin embryos. *BMC genomics*, *22*(1), 1-15.
* Enhancer RNA dataset: Khor, J. M., Guerrero-Santoro, J., Douglas, W., & Ettensohn, C. A.  (2021). Global patterns of enhancer activity during sea urchin  embryogenesis assessed by eRNA profiling. *Genome Research*, *31*(9), 1680-1692.
* lncRNA: Hezroni, H., Koppstein, D., Schwartz, M. G., Avrutin, A., Bartel, D. P., & Ulitsky, I. (2015). Principles of long noncoding RNA evolution  derived from direct comparison of transcriptomes in 17 species. *Cell reports*, *11*(7), 1110-1122.

For the above, some files contained overlaps between regions within the same file, I resolved these using the check_overlap_within_df.py.

Then, the files above along with the input_for_extra_regannot_notannot_3.1.csv file were the input to the [annot_reg_regions.py](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/annot_reg_regions.py) program. -> output list of positions that fell in an lncRNA region, and another list of positions that fell in any of the regulatory regions instead (ATAC, Chip, L.var, etc.).

#### Putting it together for final analysis

Chi-squared analysis was run as shown in chi-squared.py. It takes gathered_annotation.csv (output of gather_annot_info.py), lncdf.csv and reg_regions_withconfidence.csv (output of annot_reg_regions.py). 

To get the total number of lncRNA nucleotides (for expected number of hits in chi-squared): overlap processed with check_overlap_within_df.py, then gff annotation (processed with filtering_annotation_forlist.py to account for promoters) was subtracted with take_out_gff_from_lncrna.py. Note: sp4.lncRNAs_overlap_processed.csv was first converted from 3.1 to 5.0. Then, each region's length was calculated and summed, together with the number of nucleotides for lncRNA in ncbi. Final number (for paper lncRNA): 17,703,544, if promoter 5000: 15,343,573.

To get the total number of enhancer nucleotides (for expected number of hits in chi-squared): overlapped regions from all different resources (ATAC, Chip, etc) into 1 file using the get_enhancer_overlaps.py code (regulatory_regions_lit_review folder, output overlapped_enhancers.csv). Then, I subtracted regions in lncRNA (additional file) and gff (processed as above). Downloaded gff file for 3.1 version from NCBI, using this after processing (like above for gff 5.0). take_out_lnc_from_enhancer.py for subtracting lncRNA, then take_out_gff_from_enhancer.py for subtracting gff. Final number: 106,848,95. Alternatively, I also ran everything as above but instead of using all ATAC and DNAseq regions under different conditions (6 files in total), I only used ATAC-DNAseq overlaps (a file available in the supplementary materials). New final number: 5,185,974. Of course, I also reran annot_reg_regions.py using this alternative method, so only checking hit in ATAC-DNAseq overlap. -> this way 1.42 times more enhancer than expected instead of 1.2. If promoter 5000: 4,567,677

To get the total number of promoter nucleotides (for expected number of hits in chi-squared): calculated total number of nucleotides in gff before and after adding (or subtracting depending on the strand) 2000 bp. It was calculated this way so that overlapping is accounted for. Code is here: get_num_promoter_nuc.py. Final number: 51,056,385  (after bug fix with promoters running off chromosomes). 105,932,665 (after bug fix) if promoter region is 5000, using same code as for 2000.

Chi-squared was run in 3 different ways: one where hits in overlapping regions was counted multiple times, one where hits in overlapping regions were discarded, and one where hits in overlapping regions was counted once, as one of the regions selected at random. Results of the 3 methods were extremely similar as there were not many hits in overlapping regions.

## Gene expression



## Protein-protein interaction network



## GO enrichment

**Uniprot has GO terms** associated to each gene in the urchin genome, here is the link to retrieve that info: https://www.ebi.ac.uk/QuickGO/annotations?taxonId=7668&taxonUsage=exact. 

Used this code to transform uniprot - GO output file into the mapping file topGO expects:

```bash
awk -F "\t" '{print $2"\t"$5}' QuickGO-annotations-1642716310981-20220120.tsv > temp_mapping
sed '$!N; /^\(.*\)\n\1$/!P; D' temp_mapping > temp2_mapping # It deletes duplicate, consecutive lines from a file
awk 'BEGIN{FS="\t"} {for(i=2; i<=NF; i++) { if (!a[$1]) a[$1]=$1FS$i ;else a[$1]=a[$1]","$i};if ($1 != old) b[j++] = a[old];old=$1 } END{for (i=0; i<j; i++) print b[i] }' temp2_mapping > GO_mapping_topGO #it collapses repeated lines into 1, comma separated, output file is in this git repo
```

Now that I have all genes with associated GO terms, I need to map **LOC IDs into UniprotIDs** to get my list of interesting genes I can use in the GO enrichment together with the gene to GO mapping file above: 

To get the list of LOC: [get_loc_lists.py](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/get_loc_lists.py), takes gathered_annotation.csv, returns 3 files: LOC for promoters, LOC for non-promoters, LOC for all together. 

I can convert LOC to UniprotID using this tool: https://www.uniprot.org/uploadlists/. select From Ensemble Genomes To uniprot.

Code [LOC2Uniprot.ipynb](https://github.com/PespeniLab/urchin_local_adapt_WGS/blob/main/code/LOC2Uniprot.ipynb) to resolve cases where one LOC mapped to multiple UniprotIDs. Since most UniprotIDs sharing a single LOC mapped to the same GO terms, I just ended up selecting the the UniprotID with the most GO terms. -> list of uniprotIDs of interest, uniprotID_all_locs.txt

Using the the two key files above, I run topgo using the following code:

```R
library(topGO)

geneID2GO <- readMappings("GO_mapping_topGO") # uniprot to GO mapping
geneNames <- names(geneID2GO)

myInterestingGenes <- read.csv("uniprotIDs_prom_locs.txt", header = FALSE) # list of interesting genes, output of LOC to uniprot mapping
intgenes <- myInterestingGenes[, "V1"]
geneList <- factor(as.integer(geneNames %in% intgenes)) # mask of 0 and 1 if geneName is interesting
names(geneList) <- geneNames # geneList but annotated with the gene names

GOdata <- new("topGOdata", 
              ontology = "BP", # ontology of interest (BP, MF or CC)
              allGenes = geneList,
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher") 

allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10) # top 10 enriched terms

#showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 10, useInfo = 'all') 

```



___HERE___

### SPU for supplementary data analysis

Converted LOC (the ones in the list above, used for GO) to SPU using the GenePageGeneralInfo_AllGenes.txt file available on the echinobase. Code to make a mapping between LOC and SPU based on above file: [make_dic.py](https://github.com/Cpetak/urchin_adaptation/blob/main/code/make_dic.py) -> SPU_LOC.json Then I used this dictionary to map my LOCs to SPUs ([get_spu_list.py](https://github.com/Cpetak/urchin_adaptation/blob/main/code/get_spu_list.py)). In the case where muliple SPUs are associated to a LOC -> I kept all alternative SPUs as here I am just looking if I found a pos in a gene that is in any of the lists below. 

#### Finding interesting genes:

- Biomineralisation genes: from https://pubmed.ncbi.nlm.nih.gov/28141889/, mel_biomin.csv
- Differentially expressed genes:
  - South vs North after common garden conditions, https://onlinelibrary.wiley.com/doi/full/10.1111/evo.12036, N_S_diff_expressed.csv
  - Expression of genes of different populations in response to low pH, https://pubmed.ncbi.nlm.nih.gov/28141889
    - Genes differentially expressed between S. purpuratus populations following one day of exposure to low pH seawater, genes up-regulated in populations most frequently exposed to ph <7.8 & down-regulated in populations less frequently exposed to ph <7.8, (DE_1.csv), genes down-regulated in populations most frequently exposed to  ph <7.8 & up-regulated in populations less frequently exposed to ph <7.8 (DE_2.csv)
    - Genes differentially expressed between S. purpuratus populations following seven days of exposure to low pH seawater, genes up-regulated in populations most frequently exposed to ph  <7.8 & down-regulated in populations less frequently exposed to ph <7.8 (DE_3.csv), genes down-regulated in populations most frequently exposed to  ph <7.8 & up-regulated in populations less frequently exposed to ph  <7.8, (DE_4)
  - Review paper gathering all genes differentially expressed in low pH, https://pubmed.ncbi.nlm.nih.gov/25070868/, all_diff_expressed.csv
- Genes shown to be related to ph (i.e. SNPs correlated to pH conditions of pops), https://academic.oup.com/icb/article/53/5/857/733987, related_to_ph_snp.csv
- Artificial selection low vs normal pH, 1 vs 7 days, allele frequencies that changed, https://www.pnas.org/content/110/17/6937, 1_7_days_alleles.csv
- Genes shown to be related to ph (i.e. SNPs correlated to pH conditions during artificial selection), https://www.biorxiv.org/content/10.1101/422261v1, sig_pH75_onlySPU.txt and sig_pH80_onlySPU.txt
- mRNA restricted to the PMC lineage, https://pubmed.ncbi.nlm.nih.gov/22190640/, PMC_enriched.csv
- mRNA enriched in the PMC lineage, https://pubmed.ncbi.nlm.nih.gov/22190640/, PMC_restricted.csv
- Important genes (in GRN or master TF): TODO

Note: if I want to look for specifically TFs or things involved in ion homeostasis -> search based on GO term



## Dendrogram for biomineralisation genes

