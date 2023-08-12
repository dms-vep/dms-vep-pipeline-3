# Input data
This subdirectory contains input data used by the pipeline.

## PacBio full-length variant sequencing to link barcodes

[PacBio_amplicon.gb](PacBio_amplicon.gb): Genbank file having features to parse with [alignparse](https://jbloomlab.github.io/alignparse/). Must have *gene* (the gene of interest) and *barcode* features.

[PacBio_feature_parse_specs.yaml](PacBio_feature_parse_specs.yaml): How to parse the PacBio amplicon using [alignparse](https://jbloomlab.github.io/alignparse/).

[PacBio_runs.csv](PacBio_runs.csv): List of PacBio CCS FASTQs used to link barcodes to variants.
It must have the following columns:

 - `library`: name of the library
 - `run`: name of the sequencing run, must be unique
 - `fastq`: FASTQ file from running CCS

## Site numbering
[site_numbering_map.csv](site_numbering_map.csv): Maps sequential 1, 2, ... numbering of the gene to a "reference" numbering scheme that represents the standard naming of sites for this gene.
Also assigns each site to a region (domain) of the protein.
So must have columns *sequential_site*, *reference_site*, and *region*.

## Mutation-type classification
[data/mutation_design_classification.csv](data/mutation_design_classification.csv) classifies mutations into the different categories of designed mutations.
Should have columns *sequential_site*, *amino_acid*, and *mutation_type*.

## Neutralization standard barcodes
[neutralization_standard_barcodes.csv](neutralization_standard_barcodes.csv) barcodes for the neutralization standards.
Must have columns *barcode* and *name*, giving the barcode and name of this neutralization standard set.

## Barcode runs
[barcode_runs.csv](barcode_runs.csv) must contain the following columns (you can optionally include more):

 - `sample`: sample name, must be unique among barcode runs. Sample name must begin with `<library>-<YYMMDD>` where `<library>` is the library and `<YYMMDD>` is the date. It is recommended (but not enforced) that the full format be `<library>-<YYMMDD>-<description>-<replicate>` where `<description>` is a string description with underscores but no dashes, and `<replicate>` is a number.
 - `library`: name of library, must match a library in the barcode-variant table
 - `date`: date of sequencing, specified in a format parseable to a date by `pandas`.
 - `fastq_R1`: path to one more FASTQ R1 sequencing files, multiple files should be semicolon-delimited

## Configuration for analyzing functional effects of mutations
[func_effects_config.yml](func_effects_config.yml) has the configuration for analyzing functional effects of mutations.
The format is explained within the file.

## Configuration for analyzing antibody escape and receptor affinity
[antibody_escape_config.yml](antibody_escape_config.yml) has the configuration for analyzing effects of mutations on escape from antibodies or sera.
This same configuration can also be used to analyze escape from soluble receptor that inhibits entry to estimate how mutations affect receptor affinity.
The configuration for such receptor-inhibition experiments are also specified in this configuration as conceptually the experiments and analysis are basically the same, with just how the results are plotted and described different.
The format is explained within the file.
