#!/usr/bin/env python3
import csv
import logging
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from timeit import default_timer as timer


N_DEBUGSTEPS = 50000


logging.basicConfig(filename=snakemake.log[0],
                    format='%(asctime)s %(message)s',
                    datefmt='%d.%m.%Y %H:%M:%S | ',
                    level=logging.DEBUG)
logging.info('Started')

threshold = np.float64(snakemake.config['threshold'])
fisher_alternative = snakemake.config['fisher_alternative']
samples = pd.read_table(snakemake.input.samples_fn, index_col=False, dtype=str)
with open(snakemake.input.mgx_fn, newline='') as mgx_file:
    mgx_reader = csv.reader(mgx_file, delimiter='\t')
    next(mgx_reader)  # skip the header
    mgx = [(row[0], row[1], set(row[2:]))
           for row in mgx_reader]


start_time = timer()
with open(snakemake.output.assn_fn, 'w', newline='') as assn_file:
    assn_writer = csv.writer(assn_file, delimiter='\t')
    assn_writer.writerow([
        'feature1',
        'feature2',
        'annotation1',
        'annotation2',
        'p-value',
        'samples1',
        'samples2',
        'samples_shared'])
    with open(snakemake.input.mbx_fn, newline='') as mbx_file:
        mbx_reader = csv.reader(mbx_file, delimiter='\t')
        next(mbx_reader)  # skip the header
        steps = 0
        pairs = 0
        for row in mbx_reader:
            mb = row[0], row[1], set(row[2:])
            for mg in mgx:
                a = len(set.intersection(mb[2], mg[2]))
                b = len(mg[2]) - a
                c = len(mb[2]) - a
                d = len(samples) - a - b - c
                _, p = fisher_exact([[a, b], [c, d]],
                                    alternative=fisher_alternative)
                if p < threshold:
                    assn_writer.writerow([
                        mb[0],
                        mg[0],
                        mb[1],
                        mg[1],
                        str(p),
                        str(len(mb[2])),
                        str(len(mg[2])),
                        str(a)])
                    pairs += 1

                steps += 1
                if steps % N_DEBUGSTEPS == 0:
                    logging.debug(('{} p-values calculated, ' +
                                   '{} associations reported, ' +
                                   '{:.3f} seconds elapsed').format(
                                      steps,
                                      pairs,
                                      timer() - start_time))

logging.info(('Finished calculating {} p-values ' +
              '({} reported) in {:.3f} seconds').format(
    steps,
    pairs,
    timer() - start_time))
