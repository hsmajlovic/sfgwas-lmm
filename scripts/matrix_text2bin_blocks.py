#!/usr/bin/env python3

# converts a genotype matrix file provided as a text file
# into a binary file with 1 byte per element
# assumes elements are integers (e.g. 0,1,2,-1)
# and uses int8 for output

import math
import os
import sys

import numpy as np

# Assumes pos_qc.txt includes chromosome ID in the first column
# Variants should be sorted in the order of chromosomes (each chrom forms a contiguous block)

input_dir = sys.argv[
    1
]  # input_dir + "/partyx/combined.bin"/"count.txt"/"cov.txt"/"pheno.txt"/"snp_pos.txt" should exist
nparties = int(sys.argv[2])
counts = [
    [
        int(line)
        for line in open(os.path.join(input_dir, "party" + str(i + 1) + "/count.txt"))
    ]
    for i in range(nparties)
]
nrows = [p[0] for p in counts]
ncols = [p[1] for p in counts]
# nrow = counts[0]
# ncol = counts[1]
nfolds = int(sys.argv[3])
ncols_per_block = int(sys.argv[4])  # 8192
output_dir = sys.argv[5]

# find chrids
chr_ids = [
    line.split()[0] for line in open(os.path.join(input_dir, "party1/snp_pos.txt"))
]
ind = 0
prev = None
colchrinds = []
for ch in chr_ids:
    if prev != ch:
        colchrinds.append(ind)
    prev = ch
    ind += 1
colchrinds.append(ind)
numchr = len(colchrinds) - 1

# create paths
# ncolblock = math.ceil(float(ncol) / float(ncols_per_block))
ncolblock = 0
for ch in range(numchr):
    ncolblock += int(
        math.ceil(float(colchrinds[ch + 1] - colchrinds[ch]) / float(ncols_per_block))
    )

for p in range(nparties):
    if p == 0:
        path = os.path.join(output_dir, "party" + str(p))
        if not os.path.exists(path):
            os.makedirs(path)
    path = os.path.join(output_dir, "party" + str(p + 1))
    if not os.path.exists(path):
        os.makedirs(path)


def get_block_inds(ntot, nblock):
    blockinds = [0] * (nblock + 1)
    perblock = int(ntot / nblock)
    nrem = ntot - perblock * nblock
    for i in range(nblock):
        count = perblock
        if i < nrem:
            count += 1
        blockinds[i + 1] = blockinds[i] + count
    return blockinds


def write_block_inds(outfile, blockinds):
    sizefile = open(outfile, 'wt')
    for i in range(len(blockinds) - 1):
        sizefile.write(str(blockinds[i + 1] - blockinds[i]) + "\n")
    sizefile.close()


outfile = [[] for p in range(nparties)]
for p in range(nparties):
    for k in range(nfolds):
        arr = [
            os.path.join(
                output_dir,
                "party" + str(p + 1),
                "fold" + str(k + 1) + "." + str(j) + ".bin",
            )
            for j in range(ncolblock)
        ]
        outfile[p].append(arr)
print("ncolblock", ncolblock)

foldsizes = [np.zeros(nfolds) for i in range(nparties)]
foldpartyinds = [None for i in range(nparties)]
for p in range(nparties):
    rowpartyinds = get_block_inds(nrows[p], nfolds)
    print("party", p, rowpartyinds)
    for f in range(nfolds):
        foldsizes[p][f] = rowpartyinds[f + 1] - rowpartyinds[f]
    foldpartyinds[p] = rowpartyinds
    print("foldsizes", p, foldsizes)
    print("rowpartyinds", p, rowpartyinds)

for p in range(nparties):
    if p == 0:
        txt = "\n".join([str(int(v)) for v in foldsizes[p]])
        fp = open(os.path.join(output_dir, "party" + str(p), "foldSizes.txt"), 'wt')
        fp.write(txt)
        fp.close()
    txt = "\n".join([str(int(v)) for v in foldsizes[p]])
    fp = open(os.path.join(output_dir, "party" + str(p + 1), "foldSizes.txt"), 'wt')
    fp.write(txt)
    fp.close()

# rowpartyinds = get_block_inds(nrow, nfolds)
# rowblockinds = np.zeros(nparties*nfolds+1, dtype=np.int)
# foldsizes = [np.zeros(nfolds) for i in range(nparties)]
# ind = 1
# #for p in range(nparties):
# for p in range(nfolds):
#     n2 = rowpartyinds[p+1] - rowpartyinds[p]
#     #foldinds = get_block_inds(n2, nfolds)
#     foldinds = get_block_inds(n2, nparties)
#     for f in range(nparties):
#         foldsizes[f][p] = foldinds[f+1]-foldinds[f]
#         rowblockinds[ind] = foldinds[f+1]-foldinds[f]+rowblockinds[ind-1]
#         ind += 1

# for p in range(nparties):
#     if p == 0:
#         txt = "\n".join([str(int(v)) for v in foldsizes[p]])
#         fp = open(os.path.join(output_dir, "party"+str(p), "foldSizes.txt"), 'wt')
#         fp.write(txt)
#         fp.close()
#     txt = "\n".join([str(int(v)) for v in foldsizes[p]])
#     fp = open(os.path.join(output_dir, "party"+str(p+1), "foldSizes.txt"), 'wt')
#     fp.write(txt)
#     fp.close()

# print("fold inds")
# print(len(rowpartyinds))
# print(rowpartyinds)
# print(np.array(rowpartyinds[1:])-np.array(rowpartyinds[0:-1]))

# print("ind block inds")
# print(len(rowblockinds))
# print(rowblockinds)
# print(rowblockinds[1:]-rowblockinds[0:-1])

# colblockinds = np.zeros(ncolblock+1, dtype=np.int)
# for b in range(ncolblock):
#    colblockinds[b+1] = colblockinds[b] + ncols_per_block
#    if colblockinds[b+1] > ncol:
#        colblockinds[b+1] = ncol

block2chr = []

colblockinds = np.zeros(ncolblock + 1, dtype=int)
ind = 1
for ch in range(numchr):
    n2 = colchrinds[ch + 1] - colchrinds[ch]
    nblock = int(math.ceil(float(n2) / float(ncols_per_block)))
    for b in range(nblock):
        colblockinds[ind] = colblockinds[ind - 1] + ncols_per_block
        ind += 1
        block2chr.append(str(ch))
    colblockinds[ind - 1] = colchrinds[ch + 1]

print("chr inds")
print(len(colchrinds))
print(colchrinds)
print(np.array(colchrinds[1:]) - np.array(colchrinds[0:-1]))

print("snp block inds")
print(len(colblockinds))
print(colblockinds)
print(colblockinds[1:] - colblockinds[0:-1])

# create block files
for p in range(nparties):
    if p == 0:
        write_block_inds(
            os.path.join(output_dir, "party" + str(p), "blockSizes.txt"), colblockinds
        )
        f = open(os.path.join(output_dir, "party" + str(p), "blockToChrom.txt"), 'wt')
        f.write("\n".join(block2chr) + "\n")
        f.close()

    write_block_inds(
        os.path.join(output_dir, "party" + str(p + 1), "blockSizes.txt"), colblockinds
    )
    f = open(os.path.join(output_dir, "party" + str(p + 1), "blockToChrom.txt"), 'wt')
    f.write("\n".join(block2chr) + "\n")
    f.close()

infiles = [
    open(
        os.path.join(input_dir, "party" + str(i + 1) + "/combined_filtered.bin"),
        'rb',
    )
    for i in range(nparties)
]

for p in range(nparties):
    arr = np.fromfile(infiles[p], dtype=np.int8)
    arr = np.reshape(arr, (nrows[p], ncols[p]))
    print("shape", p, arr.shape)
    for k in range(nfolds):
        for b in range(ncolblock):
            print("processing ", outfile[p][k][b])
            f = open(outfile[p][k][b], 'ab')
            arr[
                int(foldpartyinds[p][k]) : int(foldpartyinds[p][k + 1]),
                colblockinds[b] : colblockinds[b + 1],
            ].tofile(f)
            f.close()
