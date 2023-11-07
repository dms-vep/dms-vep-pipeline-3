"""Get functional scores from barcode counts."""


import sys

import alignparse.utils

import Bio.SeqIO

import dms_variants.codonvarianttable

import pandas as pd

sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Computing functional scores for selection={snakemake.wildcards.selection}\n")

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

print(f"\nWriting summary counts to {snakemake.output.count_summary}")
counts_summary_d = {
    "selection": snakemake.wildcards.selection,
    "library": library,
    "pre_selection_sample": snakemake.params.samples["pre_selection_sample"],
    "post_selection_sample": snakemake.params.samples["post_selection_sample"],
    "pre_selection_date": snakemake.params.dates["pre_selection_sample"],
    "post_selection_date": snakemake.params.dates["post_selection_sample"],
    "min_pre_selection_count": min_pre,
}
for ctype in ["pre_count", "post_count"]:
    counts_summary_d[f"{ctype}_median"] = func_scores[ctype].median()
    counts_summary_d[f"{ctype}_q1"] = func_scores[ctype].quantile(0.25)
    counts_summary_d[f"{ctype}_q3"] = func_scores[ctype].quantile(0.75)
    counts_summary_d[f"{ctype}_min"] = func_scores[ctype].min()
    counts_summary_d[f"{ctype}_max"] = func_scores[ctype].max()
(
    pd.Series(counts_summary_d)
    .to_frame()
    .transpose()
    .to_csv(snakemake.output.count_summary, index=False, float_format="%.4g")
)

func_scores = func_scores.query("pre_count >= @min_pre")

# renumber aa_substitutions into reference numbering, also keep original numbering
# as aa_subsitutions_sequential
renumber = alignparse.utils.MutationRenumber(
    number_mapping=pd.read_csv(snakemake.input.site_numbering_map),
    old_num_col="sequential_site",
    new_num_col="reference_site",
    wt_nt_col=None,
    allow_letter_suffixed_numbers=True,
)

func_scores = (
    func_scores.query("target == 'gene'")
    .rename(
        columns={
            "aa_substitutions": "aa_substitutions_sequential",
            "codon_substitutions": "codon_substitutions_sequential",
        }
    )
    .assign(
        aa_substitutions=lambda x: x["aa_substitutions_sequential"].apply(
            renumber.renumber_muts,
            allow_gaps=True,
            allow_stop=True,
        )
    )
)

print(f"\nWriting functional scores to {snakemake.output.func_scores}")

(
    func_scores.sort_values(["func_score", "barcode"], ascending=[False, True])[
        [
            "func_score",
            "func_score_var",
            "barcode",
            "aa_substitutions",
            "n_aa_substitutions",
            "n_codon_substitutions",
        ]
    ].to_csv(snakemake.output.func_scores, float_format="%.4g", index=False)
)
