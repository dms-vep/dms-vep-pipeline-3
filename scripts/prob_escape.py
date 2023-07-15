"""Get probabilities (fraction) of antibody escape from variant counts."""


import sys

import Bio.SeqIO

import alignparse.utils

import dms_variants.codonvarianttable

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

geneseq = str(Bio.SeqIO.read(snakemake.input.gene_sequence_codon, "fasta").seq)

variants = dms_variants.codonvarianttable.CodonVariantTable(
    barcode_variant_file=snakemake.input.codon_variants,
    geneseq=geneseq,
    allowgaps=True,
    substitutions_are_codon=True,
    primary_target="gene",
    substitutions_col="codon_substitutions",
)

library = set(snakemake.params.libraries.values())
if len(library) != 1:
    raise ValueError("Samples not from 1 library: {snakemake.params.libraries}")
library = list(library)[0]

dates = set(snakemake.params.dates.values())
if len(dates) != 1:
    raise ValueError("Samples not from 1 date: {snakemake.params.dates}")

variant_counts = pd.concat(
    [
        pd.read_csv(snakemake.input.no_antibody_sample).assign(sample="no_antibody"),
        pd.read_csv(snakemake.input.antibody_sample).assign(sample="antibody"),
    ]
).assign(library=library)

variants.add_sample_counts_df(variant_counts)

prob_escape, neut_standard_fracs, _ = variants.prob_escape(
    selections_df=pd.DataFrame(
        {
            "library": [library],
            "antibody_sample": "antibody",
            "no-antibody_sample": "no_antibody",
        }
    ),
    neut_standard_target=snakemake.params.neut_standard,
    primary_target_only=True,
    min_neut_standard_frac=0,
    min_neut_standard_count=0,
)

# renumber aa_substitutions in reference numbering
renumber = alignparse.utils.MutationRenumber(
    number_mapping=pd.read_csv(snakemake.input.site_numbering_map),
    old_num_col="sequential_site",
    new_num_col="reference_site",
    wt_nt_col=None,
    allow_letter_suffixed_numbers=True,
)
prob_escape = prob_escape.rename(
    columns={"aa_substitutions": "aa_substitutions_sequential"}
).assign(
    aa_substitutions=lambda x: x["aa_substitutions_sequential"].apply(
        renumber.renumber_muts, allow_gaps=True, allow_stop=True
    )
)

# write output
prob_escape[
    [
        "barcode",
        "aa_substitutions",
        "prob_escape",
        "prob_escape_uncensored",
        "antibody_count",
        "no-antibody_count",
    ]
].to_csv(snakemake.output.prob_escape, index=False, float_format="%.4f", na_rep="nan")

neut_standard_fracs[
    ["antibody_count", "antibody_frac", "no-antibody_count", "no-antibody_frac"]
].to_csv(snakemake.output.neut_standard_fracs, index=False, float_format="%.4g")
