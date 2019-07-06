configfile: 'config.yaml'
localrules: all


rule all:
    input:
        'output/associations.tsv'


rule decoy_database:
    input:
        mbx_fn='input/mbx.feat.bin_sparse',
        samples_fn='input/samples.tsv'
    output:
        mbx_decoy_fn='input/mbx_decoy.feat.bin_sparse'
    script:
        'scripts/decoy_database.py'



if config['associations_parts'] == 1:
    # SERIAL EXECUTION
    rule associations_serial:
        input:
            mbx_fn='input/mbx.feat.bin_sparse',
            mgx_fn='input/mgx.feat.bin_sparse',
            samples_fn='input/samples.tsv'
        output:
            assn_fn=protected('output/associations.tsv')
        log:
            'log/associations.log'
        script:
            'scripts/associations.py'


    rule associations_decoy_serial:
        input:
            mbx_fn='input/mbx_decoy.feat.bin_sparse',
            mgx_fn='input/mgx.feat.bin_sparse',
            samples_fn='input/samples.tsv'
        output:
            assn_fn=protected('output/associations_decoy.tsv')
        log:
            'log/associations_decoy.log'
        script:
            'scripts/associations.py'
else:
    # PARALLEL EXECUTION
    PARTS = list(range(int(config['associations_parts'])))


    rule split_input:
        input:
            mbx_fn='input/mbx.feat.bin_sparse',
            mbx_decoy_fn='input/mbx_decoy.feat.bin_sparse',
            mgx_fn='input/mgx.feat.bin_sparse'
        output:
            mbx_parts=temp(expand('tmp/mbx.part{part}.feat.bin_sparse', part=PARTS)),
            mbx_decoy_parts=temp(expand('tmp/mbx_decoy.part{part}.feat.bin_sparse', part=PARTS)),
            mgx_parts=temp(expand('tmp/mgx.part{part}.feat.bin_sparse', part=PARTS))
        script:
            'scripts/split_input.py'


    rule associations_part:
        input:
            mbx_fn='tmp/mbx.part{part}.feat.bin_sparse',
            mgx_fn='tmp/mgx.part{part}.feat.bin_sparse',
            samples_fn='input/samples.tsv'
        output:
            assn_fn=temp('tmp/associations.part{part}.tsv')
        log:
            'log/associations.part{part}.log'
        script:
            'scripts/associations.py'


    rule associations_multiprocess:
        input:
            assn_parts=expand('tmp/associations.part{part}.tsv', part=PARTS)
        output:
            assn_fn=protected('output/associations.tsv')
        shell:
            '''
            tail_start=1
            for part in {input.assn_parts}
            do
                tail -n+$tail_start $part >>{output.assn_fn}
                tail_start=2
            done
            '''


    rule associations_decoy_part:
        input:
            mbx_fn='tmp/mbx_decoy.part{part}.feat.bin_sparse',
            mgx_fn='tmp/mgx.part{part}.feat.bin_sparse',
            samples_fn='input/samples.tsv'
        output:
            assn_fn=temp('tmp/associations_decoy.part{part}.tsv')
        log:
            'log/associations_decoy.part{part}.log'
        script:
            'scripts/associations.py'


    rule associations_decoy_multiprocess:
        input:
            assn_parts=expand('tmp/associations_decoy.part{part}.tsv', part=PARTS)
        output:
            assn_fn=protected('output/associations_decoy.tsv')
        shell:
            '''
            tail_start=1
            for part in {input.assn_parts}
            do
                tail -n+$tail_start $part >>{output.assn_fn}
                tail_start=2
            done
            '''

# Postprocessing


rule fdr:
    input:
        mbx_fn='input/mbx.feat.bin_sparse',
        mgx_fn='input/mgx.feat.bin_sparse',
        assn_fn='output/associations.tsv',
        assn_decoy_fn='output/associations_decoy.tsv'
    output:
        fdr_fn='output/fdr.tsv'
    script:
        'scripts/fdr.py'


rule tree_biclustering:
    input:
        assn_fn='output/associations.tsv',
        tree_fn='input/mgx.tree'
    output:
        tree_biclust_fn='output/mgx.biclust.tree',
        biclust_fn='output/biclusters.comp',
        assn_biclust_fn='output/associations.biclust.tsv'
    script:
        'scripts/tree_biclustering.py'
