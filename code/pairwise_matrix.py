import json
import pickle
from time import time
import matplotlib.pyplot as plt
import numba
import numpy as np

#from utils import hash_fun

# In this python script, we calculate pairwise differences between individuals and plot it as a matrix

def transpose(l):
    return list(zip(*l))

#### Load data

# (137, 991430)
# 0 : reference
# 1 : one reference, one alternative
# 2 : two alternative
# 9 : missing data
# the data has been filtered by a great wise biologist
genotypes = np.load("data/all_snp.npy")

# { chromosome : length }
with open("data/chr_lengths", "r") as f:
    chr_lengths = json.loads(f.read())

# Chromosome, position information for each of the columns in all_snp (137 x 900,000)
# ( chromosome , position ) x 991430 entries
with open("data/coords", "r") as f:
    coords = json.loads(f.read())

#####################

### Quick check on missing data
missing_per_col = (genotypes == 9).sum(axis=0)
plt.plot(np.sort(missing_per_col))
plt.savefig("fig1") # y-axis: how many chromosomes have missing data at pos (out of 280-6), x-axis: positions (max ~1e6)
# shows that at most there are around 18 missing, most of them miss around 6

### Quick check on minor allele frequency
plt.figure(figsize=(13, 10))
plt.plot(np.sort(np.where(genotypes == 9, 0, genotypes).sum(axis=0)) / (len(genotypes) * 2),",")
plt.title("Minor allele frequency over the positions")
plt.xlabel("Loci")
plt.ylabel("Frequency")
plt.savefig("fig2")

### Dividing into pops
pop_sizes = [20, 20, 20, 20, 20, 19, 18] # order is South to North, 3 outliers taken out
assert sum(pop_sizes) == 137
starts = np.cumsum([0] + pop_sizes)
ends = np.cumsum(pop_sizes)
pops = [genotypes[start:end] for start, end in zip(starts, ends)]

#####################

### Calculate pairwise distances and save in file

@numba.njit("f8(i2[:],i2[:])")
def distance(ind1, ind2):
    N = len(ind1)
    differences = 0
    tot = 0
    for i in range(N):
        allele1, allele2 = ind1[i], ind2[i]
        if allele1 != 9 and allele2 != 9:
            # abs(0 - 0) = 0
            # abs(1 - 1) = 0
            # abs(2 - 2) = 0
            #
            # abs(0 - 1) = 1
            # abs(1 - 0) = 1
            #
            # abs(1 - 2) = 1
            # abs(2 - 1) = 1
            #
            # abs(2 - 0) = 2
            # abs(0 - 2) = 2
            differences += abs(allele1 - allele2)
            tot += 1
    distance = differences / (2 * tot)  # because now the max is 2 * len
    return distance


@numba.njit(parallel=True)
def compute_distances(genotypes, distance):
    N, K = genotypes.shape
    distances = np.zeros((N, N), dtype=np.float64)
    for i in numba.prange(N):
        for j in range(i):
            if i != j:
                d = distance(genotypes[i], genotypes[j])
                distances[i, j] = d
                distances[j, i] = d
    return distances

#distances = compute_distances(genotypes, distance)
#np.save(f"data/{hash_fun(distance)}", distances)

distances = np.load(f"data/distance_7695cdaa9ce47ecd.npy")

#####################

### Plot pairwise distances

# check that only the diagonal is zero, then make it equal to the mean for plotting purposes
diag_indices = np.diag_indices(distances.shape[0])
assert np.all(diag_indices == np.argwhere(distances == 0.0).T)
mean = distances[np.nonzero(distances)].mean()
distances[diag_indices] = mean

fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(distances, interpolation="nearest")
for end in ends:
    ax.axvline(end - 1 + 0.5, lw=0.5, alpha=0.5, color="red")
    ax.axhline(end - 1 + 0.5, lw=0.5, alpha=0.5, color="red")
names = ["SAN", "LOM", "TER", "BOD", "KIB", "CAP", "FOG"]
pos = [end - 10 for end in ends]
ax.set_xticks(pos, names)
ax.set_yticks(pos, names)
ax.xaxis.tick_top()
fig.colorbar(im, ax=ax)
plt.savefig("all_pairwise.png")

### Plot averages of populations (e.g. SAN average distance from LOM, etc)

borders = list(zip(starts, ends))

def get_block(distances, pop1, pop2):
    return distances[pop1[0] : pop1[1], pop2[0] : pop2[1]]

M = len(borders)
pop_distances = np.zeros((M, M))
for i in range(M):
    for j in range(M):
        pop1 = borders[i]
        pop2 = borders[j]
        block = get_block(distances, pop1, pop2)
        pop_distances[i, j] = block.mean()

fig, ax = plt.subplots(figsize=(10, 10))
im = ax.imshow(pop_distances)
names = ["SAN", "LOM", "TER", "BOD", "KIB", "CAP", "FOG"]
pos = np.arange(len(names))
ax.set_xticks(pos, names)
ax.set_yticks(pos, names)
fig.colorbar(im, ax=ax)
plt.savefig("averaged_pairwise.png")
