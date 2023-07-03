"""Count variants from Illumina barcodes."""


import os
import sys

import dms_variants.illuminabarcodeparser

import pandas as pd


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

# get variables from Snakemake
fastq_R1 = snakemake.input.fastq_R1
variants_csv = snakemake.input.variants
counts_csv = snakemake.output.counts
counts_invalid_csv = snakemake.output.invalid
fates_csv = snakemake.output.fates
sample = snakemake.wildcards.sample
library = snakemake.params.library
parser_params = snakemake.params.parser_params

# get valid barcodes
variants = pd.read_csv(variants_csv)

if library not in variants["library"].unique():
    raise ValueError(f"{library=} not in {variants['library'].unique()=}")

valid_barcodes = variants.query("library == @library")["barcode"].tolist()
assert len(valid_barcodes) == len(set(valid_barcodes))

print(f"There are {len(valid_barcodes)} valid barcodes for {library=}")
if len(valid_barcodes) < 1:
    raise ValueError(f"no valid barcodes for {library=}")

# get barcode length
bclen = len(valid_barcodes[0])
if not all(bclen == len(bc) for bc in valid_barcodes):
    raise ValueError("not all barcodes of same length")

valid_barcodes = set(valid_barcodes)

print(f"Parsing barcodes from {fastq_R1}")
parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
    bclen=bclen,
    **parser_params,
)

for f in fastq_R1:
    if not os.path.isfile(f):
        raise ValueError(f"Cannot find FASTQ {f}")

counts, fates = parser.parse(fastq_R1)

counts["valid"] = counts["barcode"].isin(valid_barcodes)

print("Counts of valid and invalid barcodes:")
print(
    counts.groupby("valid").aggregate(
        n_barcodes=pd.NamedAgg("barcode", "nunique"),
        n_counts=pd.NamedAgg("count", "sum"),
    )
)

print(f"Writing valid barcode counts to {counts_csv}")
counts_valid = counts.query("valid").drop(columns="valid")
missing_valid_barcodes = sorted(list(valid_barcodes - set(counts_valid["barcode"])))
# add zero counts for missing valid barcodes
counts_valid = counts_valid = pd.concat(
    [
        counts_valid[["barcode", "count"]],
        pd.DataFrame({"barcode": missing_valid_barcodes, "count": 0}),
    ]
)
counts_valid.to_csv(counts_csv, index=False)

print(f"Writing invalid barcode counts to {counts_invalid_csv}")
counts_invalid = (
    counts.query("not valid")
    .drop(columns="valid")
    .assign(library=library, sample=sample)
)
counts_invalid[["barcode", "count"]].to_csv(counts_invalid_csv, index=False)

print(f"Writing barcode fates to {fates_csv}")
fates = pd.concat(
    [
        fates.query("fate not in ['valid barcode', 'invalid barcode']"),
        pd.DataFrame(
            {
                "fate": ["valid barcode", "invalid barcode"],
                "count": [
                    counts_valid["count"].sum(),
                    counts_invalid["count"].sum(),
                ],
            }
        ),
    ]
).assign(library=library, sample=sample)
fates[["fate", "count"]].to_csv(fates_csv, index=False)
