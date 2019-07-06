import csv
import pandas as pd
import numpy as np


# algorithm by Floyd, see doi:10.1145/30401.315746
def random_subset(n, k):
    s = set()
    for i in range(n - k, n):
        a = np.random.randint(i + 1)
        s.add(a if a not in s else i)
    return np.fromiter(s, dtype=np.int)


if 'seed' in snakemake.config:
    np.random.seed(snakemake.config['seed'])

samples = pd.read_table(snakemake.input.samples_fn,
                        index_col=False,
                        dtype=str).name.values

with open(snakemake.input.mbx_fn, newline='') as mbx_file, \
        open(snakemake.output.mbx_decoy_fn, 'w', newline='') as mbx_decoy_file:
    mbx_reader = csv.reader(mbx_file, delimiter='\t')
    mbx_decoy_writer = csv.writer(mbx_decoy_file, delimiter='\t')
    mbx_decoy_writer.writerow(next(mbx_reader))
    for row in mbx_reader:
        p = random_subset(len(samples), len(row) - 2)
        np.random.shuffle(p)
        mbx_decoy_writer.writerow([
            row[0],
            'MBX_DECOY' + row[1][3:],
            *[samples[i] for i in p]])
