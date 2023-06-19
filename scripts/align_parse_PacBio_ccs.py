"""Implements ``snakemake`` rule `align_parse_PacBio_ccs`."""

import os
import sys

import alignparse.minimap2
import alignparse.targets

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

os.makedirs(snakemake.output.outdir, exist_ok=True)

# do alignments and parsing
targets = alignparse.targets.Targets(
    seqsfile=snakemake.input.amplicon, feature_parse_specs=snakemake.input.specs
)
mapper = alignparse.minimap2.Mapper(alignparse.minimap2.OPTIONS_CODON_DMS)
samfile = os.path.join(snakemake.output.outdir, "alignments.sam")
targets.align(snakemake.input.fastq, samfile, mapper)
readstats, aligned, filtered = targets.parse_alignment(
    samfile,
    to_csv=True,
    csv_dir=snakemake.output.outdir,
)

# write read stats to CSV
readstats.to_csv(os.path.join(snakemake.output.outdir, "readstats.csv"), index=False)

# write CSVs giving files with aligned and filtered information
(
    pd.DataFrame.from_records(
        list(aligned.items()), columns=["target", "csv_file"]
    ).to_csv(os.path.join(snakemake.output.outdir, "aligned.csv"), index=False)
)
(
    pd.DataFrame.from_records(
        list(filtered.items()), columns=["target", "csv_file"]
    ).to_csv(os.path.join(snakemake.output.outdir, "filtered.csv"), index=False)
)
