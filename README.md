# Metagenomic-metabolomic associations workflow

## Short description
This workflow is designed to study large arrays of multiomics data by analyzing patterns of cooccurrence
of pairs of metagenomics and metabolomics features in a large amount of samples in order to find statistically
significant functional associations, e.g. a bacterium produces an antibiotic. The workflow is also meant
to assess various metabolomics feature extraction tools and/or correlation methods by generating a decoy database and
calculating False Discovery Rate (FDR). The workflow is built and run with
[Snakemake](https://snakemake.readthedocs.io/en/stable/), which allows to easily change and expand the pipeline when
necessary.

## Directories structure
* `input`: _preprocessed_ initial data
* `log`: workflow steps logs
* `output`: results
* `preprocessing`: subworkflow for converting the raw data to the unified format
* `scripts`: workflow scripts
* `tmp`: temporary workflow files

## Main steps

### Data preprocessing
Features data comes in a variety of format, depending on the source and methods used to obtain it, therefore a
preprocessing step is needed to convert this data to a unified format (described below). There are no restrictions on
the way this step is performed as long as it follows the output requirements, however, a Snakemake-based subworkflow
located in the `preprocessing` directory can be used (needs to be implemented explicitly).

As the result of this step the following files should be created in the `input` directory:

* `mbx.feat.<type>`: metabolomics features
* `mgx.feat.<type>`: metagenomics features
* `samples.tsv`: samples' names and additional information

### Decoy database creation
To be able to evaluate the quality of the feature extraction or association methods that are used, we generate a decoy
database of metabolomic features by applying a random permutation to the sample set for every feature. This decoy
database is later use to estimate the FDR.

As the result of this step the following files are created in the `input` directory:

* `mbx_decoy.feat.<type>`: decoy metabolomics features

### Pairwise associations search
For every pair of a metagenomic and a metabolomic features we compare their abundance profiles to find pairs of
associated features. For now, Fisher's exact test
([scipy.stats.fisher_exact](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.fisher_exact.html)) is
applied to binary abundance data (present or not present) and the p-value is evaluated to measure the statistical
significance of the association of these features, however, other association evaluation methods, e.g. Pearson or
Spearman correlation, may be added in the future.

Fisher's exact test can use one of the three alternative hypotheses. Depending on which `fisher_alternative` is used
different kinds of associations would be searched for:

* `greater`: positive associations (presence of feature `A` is associated with presence of feature 'B') would be found
* `less`: negative associations (presence of feature `A` is associated with absence of feature 'B') would be found
* `two-sided`: both kinds of associations would be found

Only the pairs with a p-value under the given `threshold` are reported. If `associations_parts` equals 1 then all input
pairs are processed as one job, otherwise the input is split into given amount of parts of almost equal size and each
part is processed as a separate job, all results being merged afterwards. As the result of this step the following
files are created in the `output` directory:

* `associations.tsv`: the reported pairs of associated features

### Decoy associations search
The same procedure described above is applied to the decoy database of metabolomic features to determine known false
positive associations.

The same `fisher_alternative`, `threshold` and `associations_parts` are used. As the result of this step the following
files are created in the `output` directory:

* `associations_decoy.tsv`: the reported pairs of associated features (using decoy metabolomics feature list)

### False Discovery Rate estimation
After the associations search we assess the quality of results obtained for a particular p-value threshold by
estimating the FDR using two methods:

1. Target-Decoy Approach: `FDR = associations_decoy / associations_real`, where `associations_decoy` and
`associations_real` are the amounts of associations of metagenomic features with decoy and real metabolomics features
respectively.
2. Benjamini-Hochberg method: the upper-bound estimate on the FDR is calculated as the lowest level for which all the
reported associations are considered significant when applying the
[BH step-up procedure](https://en.wikipedia.org/wiki/False_discovery_rate#Benjamini%E2%80%93Hochberg_procedure).

FDR is estimated for a list of `fdr_thresholds`. Every value in `fdr_thresholds` should be greater than or equal to the
previously used `threshold`. As the result of this step the following files are created in the `output` directory:

* `fdr.tsv`: a table of estimated FDRs for the given p-value thresholds

## Data formats

### Samples
`samples.tsv` file is a TAB-separated table of the samples data with the first line being the header and every
consecutive line describing a sample. The table has the following columns:

* `name`: name of the sample

Other columns are allowed and will be ignored by the workflow.

### Features
Every feature file has the extension `.feat.<type>`, where `<type>` denotes the data representation type in the file
and can be one of the following:

* `bin_dense`: dense representation of which samples every the features is present in
* `bin_sparse`: sparse representation of which samples the features is present in
* `real`: real-valued presence data of the features in the samples, e.g. abundance or intensity

`*.feat.bin_sparse` files have the following structure:

* first line consists of the strings `name`, `annotation` and `samples` separated by TAB characters
* every consecutive line describes one feature and consists of the name and the annotation of the feature followed by
  the list of the samples this feature is present in, all fields being separated by TAB characters
* it is required that each feature is present in at least one sample

`*.feat.bin_dense` files have the following structure:

* first line consists of the strings `name` and `annotation` followed by the list of all the samples, all fields being
  separated by TAB characters (`'\t'`)
* every consecutive line describes one feature and consists of the name and the annotation of the feature followed by
  the list of `'0'` and `'1'` characters, meaning the absence and presence of the feature in the corresponding sample
  respectively, all fields being separated by TAB characters (`'\t'`)

`*.feat.real` files have the same structure as the `*.feat.bin_dense`, but it contains real-valued data instead of
binary `'0'`/`'1'`.

`annotation` fields are either `<type>` or `<type>/<custom_annotation>`, where `<custom_annotation>`, if available,
is the any relevant annotation of the feature, e.g. the OTU name or the result of molecular database dereplication, and
`type` is one of the following:

* `MBX_DECOY`: for metabolomic decoy features
* `MBX`: for metabolomic features
* `MGX`: for metagenomic features

### Associations

`associations.tsv` and other associations files are TAB-separated tables with the first line being the header and every
consecutive line describing a reported association pair. The table has the following columns:

* `feature1`: name of the first feature of the pair
* `feature2`: name of the second feature of the pair
* `annotation1`: annotation of the first feature of the pair
* `annotation2`: annotation of the second feature of the pair
* `p-value`: p-value of the association
* `samples1`: number of the samples the first feature is present in
* `samples2`: number of the samples the second feature is present in
* `samples_shared`: number of the samples both features are present in

Other columns are allowed and will be ignored by the workflow.

### FDR list
`fdr.tsv` file is a TAB-separated table with the first line being the header and every consecutive line describing an
FDR calculated for a p-value threshold from `fdr_thresholds`. The table has the following columns:

* `FDR`: the FDR calculated using Target-Decoy Approach
* `FDR_BH`: the FDR calculated using Benjamini-Hochberg method
* `threshold`: the p-value threshold for which the FDR is calculated
* `n_decoy_assn`: the number of associations with decoy metabolomics features having the p-value under the threshold
* `n_real_assn`: the number of associations with real metabolomics features having the p-value under the threshold

Other columns are allowed and will be ignored by the workflow.

### Configuration
`config.yaml` is the configuration file written in the [YAML](http://yaml.org/) format. The following entries should be
specified:

* `fisher_alternative`: the alternative hypothesis used in Fisher's exact test
* `threshold`: the p-value threshold for the association pair to be reported
* `fdr_thresholds`: the list of p-value thresholds for which the FDR values should be calculated
* `associations_parts`: the number of jobs the associations step is split into allowing for parallel execution

## Running the workflow
To run the whole workflow perform the following steps:

1. Prepare and do the data preprocessing step using any method you choose.
2. Run `snakemake`.

Individual steps can also be executed. For more information visit
[Snakemake](https://snakemake.readthedocs.io/en/stable/) website.

## Authors
Egor Shcherbin, Hosein Mohimani

Mohimani Lab, CMU, 2018
