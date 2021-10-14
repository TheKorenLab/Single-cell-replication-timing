#! /usr/bin/env python3

import sys
import gzip
import csv
import numpy as np
from Bio.SeqIO.QualityIO import FastqGeneralIterator

barcodes = set([ x.strip()
  for x in open('10x_barcode_whitelist.txt', 'r')
    if not ('#' in x) ])
barcode_whitelist = sorted(list(barcodes))
bc_idx = { bc:idx for (idx,bc) in enumerate(barcode_whitelist) }

bc_counts = np.zeros(len(barcode_whitelist), dtype=np.int32)
bad_count = 0
for _, seq, _ in FastqGeneralIterator(sys.stdin):
    if seq in bc_idx:
        idx = bc_idx.get(seq)
        bc_counts[idx] += 1
    else:
        bad_count += 1

if bad_count / (bad_count + bc_counts.sum()) > 0.9:
    print('Error: > 90% of barcodes are invalid.')
    sys.exit(1)

high_count_set = bc_counts[bc_counts > max(bc_counts) // 10]
min_read_count = np.percentile(high_count_set, 95) // 10

valid_barcodes = set()
for bc in barcode_whitelist:
    idx = bc_idx.get(bc)
    if bc_counts[idx] > min_read_count:
        valid_barcodes.add(bc)

np.savetxt('barcode_list.csv', list(valid_barcodes), delimiter=',', fmt='%s')
