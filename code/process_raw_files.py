import json
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

# This file processes the raw data (in raw folder) and returns them in the correct format (in data folder)

# First, fix file with genotype information for each individual
# (0 = reference, 1 = one parent has alternative, 2 = both parents have the alternative allele, 9 = missing data)
# 137 rows (no PCA outliers) or 140 rows (with outliers) and around 900,000 column which are the genotpes for the genomic positions
with open("raw/StoN_137", "r") as f:
    data = f.read()
data = data.splitlines()
rows = [list(map(int, row.split(" "))) for row in data]
data = np.array(rows, dtype=np.int16)
np.save("data/all_snp", data)

# Chromosome, position information for each of the columns in all_snp (137 x 900,000)
with open("raw/chr_pos_for_cols", "r") as f:
    lines = f.read().splitlines()
    coords = []
    for line in lines:
        ch, pos = line.split(".1")
        coords.append((ch, int(pos)))
with open("data/coords", "w") as f:
    f.write(json.dumps(coords))

# Then, we process the file containing the length of each chromosome in the dataset
with open("raw/chr_lengths", "r") as f:
    lines = f.read().splitlines()
    chr_lengths = {}
    for line in lines:
        stuff, ch, whatisthis, length = line.split(" ")
        chr_lengths[ch.split(".")[0]] = int(length)
with open("data/chr_lengths", "w") as f:
    f.write(json.dumps(chr_lengths))

# Finally, we have information about chromosome and position information about genes in the genome
# Let's see which one of our ~900,000 SNPs fall into genes
annotation_df = pd.read_csv("raw/all_cleaned_annotations.csv")
genes_coords = []
for ch, pos in tqdm(coords):
    temp = annotation_df[
        (annotation_df["chr"] == ch + ".1")
        & (annotation_df["start"] <= pos)
        & (annotation_df["stop"] >= pos)
    ]  # get rows in annotation file where our loci is between start and stop
    if len(temp) == 0:
        # print("not in gene")
        genes_coords.append([])
    else:
        gene_names = [] # SNPs could fall in multiple genes at once as some genes overlap
        for i, row in temp.iterrows():
            if row.info.startswith("ID=gene"): #otherwise it is tRNA, miRNA, etc. Everything else was already cleaned from all_cleaned_annotations.csv
                gene_name = row.info.split(";")[0][8:]
                gene_names.append(gene_name)
        id = ch + "_" + str(pos)
        genes_coords.append(gene_names)
pickle.dump(genes_coords, open("data/coords_in_genes.pkl", "wb"))

# Also, let's save some information about the genes
df = pd.read_csv("raw/all_cleaned_annotations.csv")
df = df[df["type"] == "gene"]

genes_start_stop = []
genes_info = {}
for i, (_chr, start, stop, info) in df[["chr", "start", "stop", "info"]].iterrows():
    if info.startswith("ID=gene"):
        _chr = _chr.split(".1")[0]
        info = info.split(";")[0][8:]
        genes_info[info] = {
            "chr": _chr,
            "start": int(start),
            "stop": int(stop),
        }
        genes_start_stop.append((_chr, int(start), int(stop)))


with open("data/genes_info", "w") as f:
    f.write(json.dumps(genes_info))

with open("data/genes_start_stop", "w") as f:
    f.write(json.dumps(genes_start_stop))
