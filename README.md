# Pipeline for analyzing deep mutational scanning (DMS) of viral entry proteins (VEPs)
[![Release](https://img.shields.io/github/v/release/dms-vep/dms-vep-pipeline-3?logo=github)](https://github.com/dms-vep/dms-vep-pipeline-3/releases)
[![Build Status](https://github.com/dms-vep/dms-vep-pipeline-3/actions/workflows/test.yaml/badge.svg)](https://github.com/dms-vep/dms-vep-pipeline-3/actions/workflows/test.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

## Overview
This repository contains a [snakemake](https://snakemake.readthedocs.io/) pipeline for analysis of deep mutational scanning of barcoded viral entry proteins.
This pipeline was developed by the [Bloom lab](https://jbloomlab.org/).

Note that this is version 3 of the pipeline (see the [CHANGELOG](CHANGELOG.md) for details); older versions are hosted on the no-longer maintained [https://github.com/dms-vep/dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline) repo.

In order to use the pipeline, you can include this repository as a [git submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in your own repo
Your own repo containing this submodule will have the project-specific data and code.

In other words, if the repository for your specific deep mutational scanning project is called `<my_dms_repo>`, you would add [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) as a submodule to that.
The master [Snakefile](Snakefile) for the entire analysis is then located in `dms-vep-pipeline-3/Snakefile`.
The configuration file (`config.yaml`), custom rules (`custom_rules.smk`), input data (`./data/`), results (`./results/`), and HTML documentation of results (`./docs/`) are n the top-level repo.
In other words, the directory structure would look like this:

```
<my_dms_repo>
├── dms-vep-pipeline-3 [added as git submodule]
│   ├── Snakefile [Snakefile that runs entire pipeline]
│   └── ... [other contents of dms-vep-pipeline-3]
├── README.md [README for main project]
├── custom_rules.smk [any custom rules that are adding to the pipeline]
├── config.yaml [configuration for Snakefile]
├── data [subdirectory with input data]
├── results [subdirectory with results created by pipeline]
├── docs [sphinx summary of results created by pipeline]
└── <other files / subdirectories that are part of project>
```

For this to work, the top-level `config.yaml` needs specify details for your project (see more below).

## Test example
A test example use of the pipeline is in [./test_example](test_example).
In order to make the example contained in the pipeline, the organization is a bit different for this [./test_example](test_example): it is contained as a subdirectory of the pipeline whereas for actual use of the repo you will make the pipeline a submodule of `<my_dms_repo>` as described above.
Therefore, your `config.yaml` will have different values for `pipeline_path` and `docs` as indicated in the comments in [./test_example/config.yaml](test_example/config.yaml).

Despite these differences, [./test_example](test_example) provides an example of how to set up your repo.
Running the analysis in [./test_example](test_example) performs the whole analysis for the test example and creates the sphinx rendering in [./docs](docs) which can be displayed via GitHub Pages as here: [https://dms-vep.github.io/dms-vep-pipeline-3](https://dms-vep.github.io/dms-vep-pipeline-3).

## Running the pipeline
You can run the pipeline from the top-level directory using:

    snakemake -j <n_jobs> --software-deployment-method conda -s dms-vep-pipeline-3/Snakefile

Note that the `-s dms-vep-pipeline-3/Snakefile` flag specifies to use the Snakefile in the pipeline submodule.

If have already built the `dms-vep-pipeline-3` [conda](https://docs.conda.io/) environment in [environment.yml](environment.yml), you can also just activate that environment and then run the above command without the `--software-deployment-method conda` command.

Note that if the environment is updated, all rules will re-run.
If you do not want to re-run rules if the environment is updated, specify `--rerun-triggers mtime params input code`---although if you do this, you are responsible for manually deleting any output that should be changed by the environment update.

Running the pipeline will put results in `./results/` and the HTML documentation in `./docs/`.
To display the auto-built HTML documentation via GitHub Pages, set up your repo to serve documentation via GitHub pages from the `/docs` folder of the main (or master) branch [as described here](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).
The documentation will then be at `https://dms-vep.github.io/<my_dms_repo>` (assuming you are using the [https://github.com/dms-vep](https://github.com/dms-vep) organization; otherwise replace `dms-vep` with whatever account contains your repo).

## Configuring and customizing the pipeline
The main configuration for the analysis is in the `config.yaml` file that you place in your-top level repository.
An example of this configuration file is at [test_example/config.yaml](test_example/config.yaml).
Note in order for things to work properly you need to correctly set the `pipeline_path`, `docs`, `github_repo_url`, `github_blob_url`, `description`, `year`, and `authors` keys in the `config.yaml` file.

The `config.yaml` must then contain additional keys as explained in [test_example/config.yaml](test_example/config.yaml).
The pipeline can in principle perform many types of analyses, including building codon-variant tables from PacBio sequencing, counting barcodes from experiments,  analyzing the functional effects of mutations, analyzing antibody escape, etc. Each of these steps is included in a separate file (eg, [build_variants.smk](build_variants.smk), [count_variants.smk](count_variants.smk), etc) that are included in the main [Snakefile](Snakefile).
As detailed in [test_example/config.yaml](test_example/config.yaml), you do not have to perform all of these depending on how you structure the configuration.
Analyses other than the codon-variant building and barcode counting each have their own YAML configuration which you should place in your `./data/` subdirectory.
For examples, see [test_example/data/func_effects_config.yml](test_example/data/func_effects_config.yml), [test_example/data/antibody_escape_config.yml](test_example/data/antibody_escape_config.yml), etc.
Note that these files will be easier to understand if you first familiarize yourself with the YAML language including the merge (`<<`), alias (`&`), and anchor (`*`) operators.

If you want to add custom Snakemake rules, create a file called `custom_rules.smk` in your top-level repo, and [Snakefile](Snakefile) in the `dms-vep-pipeline-3` will automatically include and run it.
For an example, see [test_example/custom_rules.smk](test_example/custom_rules.smk).
The snakemake pipeline in [Snakefile](Snakefile) defines a nested dictionary named `docs` which is used to build the HTML documentation that goes into [./docs/](docs) for rendering on GitHub pages.
So if you want other files or outputs for any of your custom rules to be rendered in these docs, add them to the nested dict `docs` as illustrated in [test_example/custom_rules.smk](test_example/custom_rules.smk).

Your top-level repo also needs to have an appropriate `.gitignore` file that specifies which results to include and which not to include (we don't include them all or it would be too much to track).
A template `.gitignore` file is at [test_example/.gitignore](test_example/.gitignore).

## Re-running the pipeline without the FASTQ files
The pipeline can perform the entire analysis starting with the FASTQ files holding the results of the PacBio sequencing (used to build the barcode-variant lookup table) and the Illumina sequencing (used to count the barcodes) to the final processed results.
However, these FASTQ files are typically very large and so are **not** tracked in the GitHub repo but need to be stored somewhere else like on a computing cluster, either at the locations they were originally generated or (for secondary use) where they were downloaded from the NCB Sequence Read Archive (SRA).
The locations of those files are specifed in the `pacbio_runs` and `barcode_runs` CSV files indicated in the `config.yaml`, and will be specific to the configuration of the computing cluster for which the pipeline is being run since the files are too large to store within a GitHub repo.

However, for many re-use purposes, secondary users do not really need to re-process the FASTQ files are the barcode counting and barcode-variant lookup table construction are fairly simple, and secondary users may just be happy to use the lookup table and counts generated from prior processing of the FASTQ files.
If you are using a repo where these counts and the barcode-variant lookup table are already computed and stored in the repo, you can then just start with those and avoid having to handle the FASTQ files at all.
To do that, you set the following options in the configuration YAML (`config.yaml`) as follows:

    prebuilt_variants: results/variants/codon_variants.csv  # use codon-variant table already in repo
    prebuilt_geneseq: results/gene_sequence/codon.fasta  # use gene sequence already in repo
    ...
    use_precomputed_barcode_counts: false  # use barcode counts already in repo

Then running the repo will no longer require any FASTQ files, and will juse utilize the precomputed variants and counts from those files stored in the repo.

## Setting up the pipeline as a submodule
To add [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline-3) as a [submodule](https://git-scm.com/book/en/v2/Git-Tools-Submodules) in your repo (`<my_dms_repo>`), do as follows.

        git submodule add https://github.com/dms-vep/dms-vep-pipeline-3

This adds the file [.gitmodules](.gitmodules) and the submodule `dms-vep-pipeline-3`, which can then be committed with:

    git commit -m 'added `dms-vep-pipeline-3` as submodule'

Note that if you want a specific commit or tag of [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline), follow the [steps here](https://stackoverflow.com/a/10916398):

    cd dms-vep-pipeline-3
    git checkout <commit>

and then `cd ../` back to the top-level directory of `<my_dms_repo>`, and add and commit the updated `dms-vep-pipeline-3` submodule.

You can also make changes to the [dms-vep-pipeline-3](https://github.com/dms-vep/dms-vep-pipeline) submodule in `<my_dms_repo>` by going into that directory, making changes on a branch, and then pushing back to [dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline-3) and opening a pull request.

## Details on setup of this pipeline repository

### Other details about pipeline

### `conda` environment
The [conda](https://docs.conda.io/) environment for the pipeline is in [environment.yml](environment.yml).

### Code formatting and linting
The Python code should be formatted with [black](https://black.readthedocs.io/) by running `black .`

Comparable formatting is done for the `snakemake` file (`*.smk` files) with [snakefmt](https://github.com/snakemake/snakefmt) by running `snakefmt .`.
The overall `snakemake` pipeline is linted by going to [./test_example](test_example) and running `snakemake -s ../Snakefile --lint`.

The code and Jupyter notebooks are linted with [ruff](https://github.com/astral-sh/ruff) by running `ruff check .` and `nbqa ruff notebooks`.

### Testing of pipeline with GitHub Actions
The pipeline is tested with GitHub Actions by checking all the formatting and linting above, and then also running the pipeline on the example in [./test_example](test_example).
See [.github/workflows/test.yaml](.github/workflows/test.yaml) for details.

### Stripping of Jupyter notebook output
The repo was configured to strip output from Jupyter notebooks as described [here](http://mateos.io/blog/jupyter-notebook-in-git/) by running:

    nbstripout --install --attributes .gitattributes

### Tracking of large files with `git lfs`
The large data files for the test example in [./test_example/sequencing_data/](test_example/sequencing_data/) are tracked with [git lfs](http://arfc.github.io/manual/guides/git-lfs).
Note the first time you set up the repo, you have to run:

    git lfs install

## Configuring a custom homepage

You can configure a custom homepage for your project with the static site generator [VitePress](https://vitepress.dev/), which can allow you to make a much nicer looking homepage.
Check out this [nicely styled homepage](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS) as an example.
Detailed instructions for adding a homepage to your project are located [here](homepage/README.md).
Note if you use this option you will need to manually edit the relevant markdown files as described in the instructions, and set GitHub Pages to serve from the `gh-pages` branch rather than from the `docs` folder (which has the automatically built documentation).
