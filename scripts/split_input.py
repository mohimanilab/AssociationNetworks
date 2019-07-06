#!/usr/bin/env python3
import csv
import numpy as np
from itertools import product


def multiples(n):
    for i in range(1, n + 1):
        if n % i == 0:
            yield i, n // i


def read_features(feat_fn):
    with open(feat_fn, newline='') as feat_file:
        feat_reader = csv.reader(feat_file, delimiter='\t')
        header = next(feat_reader)
        features = list(feat_reader)
    return features, header


mbx, _ = read_features(snakemake.input.mbx_fn)
mbx_decoy, _ = read_features(snakemake.input.mbx_decoy_fn)
mgx, header = read_features(snakemake.input.mgx_fn)

n = snakemake.config['associations_parts']
min_part_size = n + 1  # find the minimal part size to make tmp files small
for a, b in multiples(n):
    if a <= len(mbx) and b <= len(mgx) and n / min(a, b) < min_part_size:
        mbx_coords = np.linspace(0, len(mbx), num=a + 1, dtype=np.int)
        mgx_coords = np.linspace(0, len(mgx), num=b + 1, dtype=np.int)
        min_part_size = n / min(a, b)
        break
else:
    mbx_coords = np.linspace(0, len(mbx), num=n + 1, dtype=np.int)
    mgx_coords = np.linspace(0, len(mgx), num=2, dtype=np.int)
assert (len(mbx_coords) - 1) * (len(mgx_coords) - 1) == n

for (i, j), mbx_part, mbx_decoy_part, mgx_part in zip(
        product(range(len(mbx_coords) - 1),
                range(len(mgx_coords) - 1)),
        snakemake.output.mbx_parts,
        snakemake.output.mbx_decoy_parts,
        snakemake.output.mgx_parts):
    with open(mgx_part, 'w', newline='') as mgx_file, \
         open(mbx_part, 'w', newline='') as mbx_file, \
         open(mbx_decoy_part, 'w', newline='') as mbx_decoy_file:
        mbx_writer = csv.writer(mbx_file, delimiter='\t')
        mbx_decoy_writer = csv.writer(mbx_decoy_file, delimiter='\t')
        mgx_writer = csv.writer(mgx_file, delimiter='\t')
        mbx_writer.writerow(header)
        mbx_decoy_writer.writerow(header)
        mgx_writer.writerow(header)

        for k in range(mbx_coords[i], mbx_coords[i + 1]):
            mbx_writer.writerow(mbx[k])
            mbx_decoy_writer.writerow(mbx_decoy[k])
        for k in range(mgx_coords[j], mgx_coords[j + 1]):
            mgx_writer.writerow(mgx[k])
