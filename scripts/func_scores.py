"""Get functional scores from barcode counts."""


import sys

import Bio.SeqIO

import dms_variants.codonvarianttable

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(
    f"Computing functional scores for condition={snakemake.wildcards.condition}, "
    f"selection={snakemake.wildcards.selection}\n"
)

if len(set(snakemake.params.libraries.values())) != 1:
    raise ValueError(f"samples not from same library: {snakemake.params.libraries}")
library = list(snakemake.params.libraries.values())[0]
print(f"These samples are for {library=}\n")

variants = dms_variants.codonvarianttable.CodonVariantTable(
    barcode_variant_file=snakemake.input.codon_variants,
    geneseq=str(Bio.SeqIO.read(snakemake.input.gene_sequence_codon, "fasta").seq),
    allowgaps=True,
    substitutions_are_codon=True,
    primary_target="gene",
    substitutions_col="codon_substitutions",
)

for sample in ["pre_selection_sample", "post_selection_sample"]:
    barcode_csv = getattr(snakemake.input, sample)
    print(f"Getting {sample} barcode counts from {barcode_csv}\n")
    variants.addSampleCounts(library, sample, pd.read_csv(barcode_csv))

func_score_params = snakemake.params.func_score_params

print(f"Using pseudocount {func_score_params['pseudocount']}")
func_scores = variants.func_scores(
    preselection="pre_selection_sample",
    pseudocount=func_score_params["pseudocount"],
    libraries=[library],
)
assert (func_scores["pre_sample"] == "pre_selection_sample").all()
assert (func_scores["post_sample"] == "post_selection_sample").all()
assert (func_scores["library"] == library).all()

for count_type in ["pre_count", "post_count"]:
    tot = func_scores[count_type].sum()
    min_wt = max(
        func_score_params["min_wt_count"], func_score_params["min_wt_frac"] * tot
    )
    print(
        f"Requiring at least {min_wt} total wildtype counts for {count_type}, "
        f"which is the greater of min_wt_count={func_score_params['min_wt_count']} and "
        f"min_wt_frac={func_score_params['min_wt_frac']} of the {tot} total counts.\n"
    )
    wt = func_scores[f"{count_type}_wt"].unique()
    assert len(wt) == 1
    wt = wt[0]
    if wt >= min_wt:
        print(f"Adequate {count_type} wildtype counts of {wt}\n")
    else:
        raise ValueError(
            f"Insufficient {count_type} wildtype counts of {wt} (need {min_wt})"
        )

min_pre = max(
    func_score_params["min_pre_selection_count"],
    func_score_params["min_pre_selection_frac"] * func_scores["pre_count"].sum(),
)
print(f"Only keeping scores for variants with at least {min_pre} pre-selection counts:")
print(
    func_scores.assign(
        sufficient_pre_selection_counts=lambda x: x["pre_count"] >= min_pre
    )
    .groupby("sufficient_pre_selection_counts")
    .aggregate(n_variants=pd.NamedAgg("aa_substitutions", "count"))
)

func_scores = func_scores.query("pre_count >= @min_pre")

print(f"\nWriting functional scores to {snakemake.output.csv}")

(
    func_scores[
        [
            "func_score",
            "func_score_var",
            "aa_substitutions",
            "n_aa_substitutions",
            "codon_substitutions",
            "n_codon_substitutions",
        ]
    ].to_csv(snakemake.output.csv, float_format="%.4g", index=False)
)
