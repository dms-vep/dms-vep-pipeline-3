"""``snakemake`` rules for calculating antibody/sera escape or receptor affinity.

Note that this analyzes data for two assay types: 'antibody_escape' and
'receptor_affinity', corresponding to escape from antibody / sera or escape
from soluble receptor used to estimate affinity. The workflow is basically the
same, so 'antibody' is used as a synonym in some of the definitions of the
receptor affinity selections.

"""


# read the config for antibody escape
with open(config["antibody_escape_config"]) as f:
    antibody_escape_config = yaml.safe_load(f)

# get configuration for any antibody escape or receptor affinity selections
assays = ["antibody_escape", "receptor_affinity"]
assay_selections = {}
for assay in assays:
    if (sel_key := assay.split("_")[0] + "_selections") in antibody_escape_config:
        assay_selections[assay] = antibody_escape_config[sel_key]
    else:
        assay_selections[assay] = {}

#  make sure all samples defined
for sel_name, sel_d in itertools.chain.from_iterable(
    assay_selections[assay].items() for assay in assays
):
    for s in [sel_d["no_antibody_sample"], *list(sel_d["antibody_samples"])]:
        if s not in sample_to_library:
            raise ValueError(f"sample {s} for {selection_name} not in barcode_runs")

# Names and values of files to add to docs
assay_docs = {assay: collections.defaultdict(dict) for assay in assays}


wildcard_constraints:
    assay="|".join(assays),


rule prob_escape:
    """Compute probability (fraction) escape for each variant."""
    input:
        no_antibody_sample=lambda wc: os.path.join(
            "results/barcode_counts",
            assay_selections[wc.assay][wc.selection]["no_antibody_sample"]
            + "_counts.csv",
        ),
        antibody_sample="results/barcode_counts/{sample}_counts.csv",
        codon_variants=config["codon_variants"],
        gene_sequence_codon=config["gene_sequence_codon"],
        site_numbering_map=config["site_numbering_map"],
    output:
        **{
            metric: "results/{assay}/by_selection/{selection}/{sample}_"
            + metric
            + ".csv"
            for metric in ["prob_escape", "neut_standard_fracs"]
        },
    params:
        neut_standard=lambda wc: assay_selections[wc.assay][wc.selection][
            "neut_standard_name"
        ],
        # script checks that dates and libraries match for all samples
        dates=lambda wc: {
            "no_antibody_sample": sample_to_date[
                assay_selections[wc.assay][wc.selection]["no_antibody_sample"]
            ],
            "antibody_sample": sample_to_date[wc.sample],
        },
        libraries=lambda wc: {
            "no_antibody_sample": sample_to_library[
                assay_selections[wc.assay][wc.selection]["no_antibody_sample"]
            ],
            "antibody_sample": sample_to_library[wc.sample],
        },
    conda:
        "environment.yml"
    log:
        "results/logs/prob_escape_{assay}_{selection}/{sample}.txt",
    script:
        "scripts/prob_escape.py"


rule fit_escape:
    """Fit a ``polyclonal`` model to a selection."""
    input:
        spatial_distances=lambda wc: (
            [
                assay_selections[wc.assay][wc.selection]["polyclonal_params"][
                    "spatial_distances"
                ]
            ]
            if assay_selections[wc.assay][wc.selection]["polyclonal_params"][
                "spatial_distances"
            ]
            else []
        ),
        plot_hide_stats=lambda wc: [
            d["csv"]
            for d in assay_selections[wc.assay][wc.selection][
                "plot_hide_stats"
            ].values()
        ],
        prob_escapes=lambda wc: [
            rules.prob_escape.output.prob_escape.format(
                assay=wc.assay, selection=wc.selection, sample=sample
            )
            for sample in assay_selections[wc.assay][wc.selection]["antibody_samples"]
        ],
        neut_standard_fracs=lambda wc: [
            rules.prob_escape.output.neut_standard_fracs.format(
                assay=wc.assay, selection=wc.selection, sample=sample
            )
            for sample in assay_selections[wc.assay][wc.selection]["antibody_samples"]
        ],
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/fit_antibody_escape.ipynb"),
    output:
        prob_escape_mean="results/{assay}/by_selection/{selection}_prob_escape_mean.csv",
        pickle="results/{assay}/by_selection/{selection}_polyclonal_model.pickle",
        nb="results/notebooks/fit_escape_{assay}_{selection}.ipynb",
    params:
        params_yaml=lambda wc, input: yaml.dump(
            {
                "params": assay_selections[wc.assay][wc.selection],
                "neut_standard_frac_csvs": list(input.neut_standard_fracs),
                "prob_escape_csvs": list(input.prob_escapes),
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/fit_escape_{assay}_{selection}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p prob_escape_mean_csv {output.prob_escape_mean} \
            -p pickle_file {output.pickle} \
            -p selection {wildcards.selection} \
            -y "{params.params_yaml}" \
            &> {log}
        """


for assay, sels in assay_selections.items():
    for sel in sels:
        for sample in sels[sel]["antibody_samples"]:
            assay_docs[assay][
                "Probability (fraction) escape for each variant in each "
                + assay.replace("_", " ")
                + " selection (CSVs)"
            ][f"{sel} {sample}"] = rules.prob_escape.output.prob_escape.format(
                assay=assay, selection=sel, sample=sample
            )
            assay_docs[assay][
                "Fits of polyclonal models to individual "
                + assay.replace("_", " ")
                + " selections"
            ][sel] = rules.fit_escape.output.nb.format(assay=assay, selection=sel)


avg_antibody_escape_config = antibody_escape_config["avg_antibody_escape"]


rule avg_antibody_escape:
    """Average antibody escape across several selections for the same antibody/serum."""
    input:
        plot_hide_stats=lambda wc: [
            d["csv"]
            for d in avg_antibody_escape_config[wc.antibody][
                "plot_hide_stats"
            ].values()
        ],
        prob_escape_means=lambda wc: [
            rules.fit_escape.output.prob_escape_mean.format(
                assay="antibody_escape", selection=sel
            )
            for sel in avg_antibody_escape_config[wc.antibody]["selections"]
        ],
        pickles=lambda wc: [
            rules.fit_escape.output.pickle.format(
                assay="antibody_escape", selection=sel
            )
            for sel in avg_antibody_escape_config[wc.antibody]["selections"]
        ],
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/avg_antibody_escape.ipynb"),
    output:
        pickle="results/antibody_escape/averages/{antibody}_polyclonal_model.pickle",
        escape_csv="results/antibody_escape/averages/{antibody}_mut_escape.csv",
        icXX_csv="results/antibody_escape/averages/{antibody}_mut_icXX.csv",
        escape_html="results/antibody_escape/averages/{antibody}_mut_escape_nolegend.html",
        icXX_html="results/antibody_escape/averages/{antibody}_mut_icXX_nolegend.html",
        nb="results/notebooks/avg_antibody_escape_{antibody}.ipynb",
    params:
        params_yaml=lambda wc: yaml.dump(
            {"params": avg_antibody_escape_config[wc.antibody]}
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/avg_antibody_escape_{antibody}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p avg_pickle_file {output.pickle} \
            -p escape_csv {output.escape_csv} \
            -p icXX_csv {output.icXX_csv} \
            -p escape_html {output.escape_html} \
            -p icXX_html {output.icXX_html} \
            -y '{params.params_yaml}' \
            &> {log}
        """


for heading, fname in [
    ("Average selections for antibody/serum", rules.avg_antibody_escape.output.nb),
    (
        "Antibody/serum mutation escape CSVs",
        rules.avg_antibody_escape.output.escape_csv,
    ),
    ("Antibody/serum mutation ICXX CSVs", rules.avg_antibody_escape.output.icXX_csv),
    (
        "Antibody/serum mutation escape plots",
        "results/antibody_escape/averages/{antibody}_mut_escape.html",
    ),
    (
        "Antibody/serum mutation ICXX plots",
        "results/antibody_escape/averages/{antibody}_mut_icXX.html",
    ),
]:
    for antibody in avg_antibody_escape_config:
        assay_docs["antibody_escape"][heading][antibody] = fname.format(
            antibody=antibody
        )

for assay, assay_doc in assay_docs.items():
    if assay_doc:
        docs[
            {
                "antibody_escape": "Antibody/serum escape",
                "receptor_affinity": "Receptor affinity",
            }[assay]
        ] = assay_doc
