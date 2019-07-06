#!/urs/bin/env python3
from collections import defaultdict, Counter
import numpy as np
import csv


# unique_mgx = defaultdict(set)
cnt_real = Counter()
max_p = defaultdict(lambda: 0)
with open(snakemake.input.assn_fn, newline='') as assn_file:
    assn_reader = csv.reader(assn_file, delimiter='\t')
    next(assn_reader)  # skip the header
    for row in assn_reader:
        for threshold in snakemake.config['fdr_thresholds']:
            # threshold is a string here
            if np.float64(row[4]) < np.float64(threshold):
                # mg_feat = row[0] if row[2].startswith('MGX') else row[1]
                # unique_mgx[threshold].add(mg_feat)
                cnt_real[threshold] += 1
                max_p[threshold] = max(max_p[threshold], np.float64(row[4]))

# unique_mgx_decoy = defaultdict(set)
cnt_decoy = Counter()
with open(snakemake.input.assn_decoy_fn, newline='') as assn_decoy_file:
    assn_decoy_reader = csv.reader(assn_decoy_file, delimiter='\t')
    next(assn_decoy_reader)  # skip the header
    for row in assn_decoy_reader:
        for threshold in snakemake.config['fdr_thresholds']:
            # threshold is a string here
            if np.float64(row[4]) < np.float64(threshold):
                # mg_feat = row[0] if row[2].startswith('MGX') else row[1]
                # unique_mgx_decoy[threshold].add(mg_feat)
                cnt_decoy[threshold] += 1

with open(snakemake.input.mbx_fn) as mbx_file, \
         open(snakemake.input.mgx_fn) as mgx_file:
    n_total_pairs = (len(mbx_file.readlines()) - 1) * \
        (len(mgx_file.readlines()) - 1)
with open(snakemake.output.fdr_fn, 'w', newline='') as fdr_file:
    fdr_writer = csv.writer(fdr_file, delimiter='\t')
    fdr_writer.writerow([
        'FDR',
        'FDR_BH',
        'threshold',
        # 'n_decoy_mgx',
        # 'n_real_mgx',
        'n_decoy_assn',
        'n_real_assn'])
    for threshold in snakemake.config['fdr_thresholds']:
        # n_decoy_mgx = len(unique_mgx_decoy[threshold])
        # n_real_mgx = len(unique_mgx[threshold])
        n_decoy_assn = cnt_decoy[threshold]
        n_real_assn = cnt_real[threshold]
        fdr_writer.writerow([
            # str(n_decoy_mgx / n_real_mgx) if n_real_mgx else 'nan',
            str(n_decoy_assn / n_real_assn) if n_real_assn else 'nan',
            str(n_total_pairs / n_real_assn * max_p[threshold])
            if n_real_assn else 'nan',
            threshold,
            # str(n_decoy_mgx),
            # str(n_real_mgx),
            str(n_decoy_assn),
            str(n_real_assn)])
