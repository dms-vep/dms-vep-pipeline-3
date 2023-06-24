"""Top-level ``snakemake`` file that runs analysis."""


import os

import flatdict

import pandas as pd


configfile: "config.yaml"


# Read the barcode runs
barcode_run_req_cols = ["sample", "library", "date", "fastq_R1"]
if config["barcode_runs"] is None:
    # just create empty data frame with required columns
    barcode_runs = pd.DataFrame(columns=barcode_run_req_cols)
else:
    barcode_runs = pd.read_csv(config["barcode_runs"], parse_dates=["date"]).assign(
        fastq_R1=lambda x: x["fastq_R1"].str.split(";")
    )
    if not set(barcode_run_req_cols).issubset(barcode_runs.columns):
        raise ValueError(f"{barcode_runs.columns=} missing {barcode_run_req_cols=}")
assert len(barcode_runs) == barcode_runs["sample"].nunique()

# `docs` is a nested dictionary used to build HTML documentation. At the bottom of the
# nesting, keys should be short titles for files and values should be the path of the
# file created by the pipeline. Files can be any of:
#  - Jupyter notebooks (extension ".ipynb") which are converted to HTML for docs
#  - HTML files (extension ".html") such as Altair plots
#  - CSV files (extension ".csv") which are linked to on the GitHub repo
# Higher levels of nesting can be keys that are headings (or subheadings) with values
# being dictionaries that gives files and their paths.
# The pipeline rules included below automatically add to this dictionary. You can also
# make further additions of course.
docs = {}


# include pipeline rules, which also add to `docs` dictionary
include: "build_variants.smk"


if len(barcode_runs) > 0:

    include: "count_variants.smk"


# add any custom rules
custom_rules = "custom_rules.smk"
if os.path.isfile(custom_rules):

    include: os.path.join(os.path.relpath(".", config["pipeline_path"]), custom_rules)


rule all:
    input:
        flatdict.FlatDict(docs).values(),  # all files specified in `docs`
