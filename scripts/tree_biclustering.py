#!/usr/bin/env python3
import numpy as np
import pandas as pd
from skbio import TreeNode

p_threshold = np.float64(snakemake.config['tree_biclustering']['p_threshold'])
clade_threshold = \
    np.float64(snakemake.config['tree_biclustering']['clade_threshold'])

assn = pd.read_table(snakemake.input.assn_fn, index_col=False,
                     header=0,
                     dtype=str)
assn['p-value'] = assn['p-value'].astype(np.float64)
assn = assn[assn['p-value'] < p_threshold]
assn.reset_index(drop=True, inplace=True)
mbx = assn[['feature1', 'annotation1']].drop_duplicates()
mgx = assn[['feature2', 'annotation2']].drop_duplicates()

tree = TreeNode.read(snakemake.input.tree_fn)
tree = tree.shear(mgx['feature2'].values)
tree.assign_ids()
tree_size = tree.count()
clade_size = np.zeros(tree_size)
for tip in tree.tips(include_self=True):
    clade_size[tip.id] += 1
    for anc in tip.ancestors():
        clade_size[anc.id] += 1

biclust = [set() for _ in range(tree_size)]
for mb in mbx['feature1'].values:
    cnt_assn = np.zeros(tree_size)
    for mg in assn['feature2'][assn['feature1'] == mb].values:
        node = tree.find(mg)
        cnt_assn[node.id] += 1
        for anc in node.ancestors():
            cnt_assn[anc.id] += 1
    taken = set()
    for node in tree.preorder(include_self=True):
        if taken.intersection(set(anc.id for anc in node.ancestors())):
            continue
        if node.is_tip() and cnt_assn[node.id] or \
                cnt_assn[node.id] >= clade_threshold * clade_size[node.id]:
            biclust[node.id].add(mb)
            taken.add(node.id)

with open(snakemake.output.biclust_fn, 'w') as biclust_file:
    for node in tree.traverse(include_self=True):
        biclust_mbx = biclust[node.id]
        if not biclust_mbx:
            continue
        biclust_mgx = set(tip.name for tip in node.tips(include_self=True))
        name = 'Bicluster_{}_{}_{}'.format(node.id,
                                           len(biclust_mgx),
                                           len(biclust_mbx))
        biclust_file.write(name + '\n')
        for mg in mgx[mgx['feature2'].isin(biclust_mgx)].values:
            biclust_file.write('\t'.join(mg) + '\n')
        for mb in mbx[mbx['feature1'].isin(biclust_mbx)].values:
            biclust_file.write('\t'.join(mb) + '\n')
        biclust_file.write('\n')

assn_biclust = np.zeros(assn.shape[0])
for i, row in assn.iterrows():
    mb, mg = row['feature1'], row['feature2']
    node = tree.find(mg)
    while mb not in biclust[node.id]:
        node = node.parent
    assn_biclust[i] = node.id
assn.assign(bicluster_id=assn_biclust).to_csv(snakemake.output.assn_biclust_fn,
                                              sep='\t',
                                              index=False)

for node in tree.non_tips(include_self=True):
    node.name = str(node.id) + ' (' + str(len(biclust[node.id])) + ')'
for node in tree.tips():
    node.name = node.name + ' (' + str(len(biclust[node.id])) + ')'
tree.write(snakemake.output.tree_biclust_fn)
