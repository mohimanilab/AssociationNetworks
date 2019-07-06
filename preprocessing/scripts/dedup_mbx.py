#!/usr/bin/env python3
import csv
import logging
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import fisher_exact
from timeit import default_timer as timer


logging.basicConfig(filename=snakemake.log[0],
                    format='%(asctime)s %(message)s',
                    datefmt='%d.%m.%Y %H:%M:%S | ',
                    level=logging.DEBUG)
logging.info('Started')

mz_list = \
    pd.read_table(snakemake.input.mz_list_fn, index_col=0)['m/z'].sort_values()
mbx = dict()
with open(snakemake.input.mbx_fn, newline='') as mbx_file:
    mbx_reader = csv.reader(mbx_file, delimiter='\t')
    next(mbx_reader)
    for row in mbx_reader:
        mbx[row[0]] = row[1][4:], set(row[2:])

start_time = timer()
graph = nx.Graph()
graph.add_nodes_from(mbx)
n_features = len(mbx)
n_samples = len(pd.read_table(snakemake.input.mbx_samples_fn,
                              index_col=False).name)
pairs = 0
for i in range(n_features):
    feat1 = mz_list.index[i]
    mz1 = mz_list[i]
    for j in range(i + 1, n_features):
        feat2 = mz_list.index[j]
        mz2 = mz_list[j]
        if mz2 - mz1 >= np.float64(snakemake.config['mz_threshold']):
            break
        pairs += 1
        a = len(set.intersection(mbx[feat1][1], mbx[feat2][1]))
        b = len(mbx[feat2][1]) - a
        c = len(mbx[feat1][1]) - a
        d = n_samples - a - b - c
        _, p = fisher_exact([[a, b], [c, d]], alternative='greater')
        if p < np.float64(snakemake.config['p_threshold']):
            graph.add_edge(feat1, feat2)
        logging.debug(('Features {} and {} tested, ' +
                       'm/z diff: {:.3f}, ' +
                       'p-value: {}, {}').format(
                          feat1,
                          feat2,
                          mz2 - mz1,
                          p,
                          'ADDED'
                          if p < np.float64(snakemake.config['p_threshold'])
                          else 'DISCARDED'))
logging.info(('Finished building the graph in {:.3f} seconds,' +
              '{} edges out of {} candidate pairs were added').format(
                  timer() - start_time,
                  graph.number_of_edges(),
                  pairs))

with open(snakemake.output.mbx_res_fn, 'w', newline='') as mbx_res_file, \
        open(snakemake.output.mbx_comp_fn, 'w') as mbx_comp_file:
    mbx_res_writer = csv.writer(mbx_res_file, delimiter='\t')
    mbx_res_writer.writerow(['feature', 'annotation', 'samples'])
    for i, comp in enumerate(nx.connected_components(graph)):
        mz_mean = mz_list[comp].mean()
        name = '_'.join(('Comp',
                         str(i),
                         str(len(comp)),
                         '{:.3f}'.format(mz_mean)))
        anno_set = set(mbx[feat][0] for feat in comp if mbx[feat][0])
        anno_policy = snakemake.config['dedup_annotations']
        if anno_policy == 'union' and len(anno_set) > 0 or \
                anno_policy == 'consensus' and len(anno_set) == 1:
                    anno = 'MBX/' + '###'.join(anno_set)
        else:
            anno = 'MBX'
        mbx_res_writer.writerow([name, anno, *set.union(*[mbx[feat][1]
                                                          for feat in comp])])

        mbx_comp_file.write(name + '\n')
        for feat in comp:
            mbx_comp_file.write(feat + '\t' + mbx[feat][0] + '\n')
        mbx_comp_file.write('\n')
logging.info('Wrote the deduplicated features to the output')
