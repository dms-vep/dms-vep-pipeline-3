"""``snakemake`` rules for escape from some agent in assays (eg, antibody escape).

Note that these rules can analyze arbitrary assays. Throughout the workflow, the term
'antibody' is used as a synonym for the agent that is being used to select for escape.

"""


# read the config for antibody escape
with open(config["antibody_escape_config"]) as f:
    antibody_escape_config = yaml.safe_load(f)

# get configuration for any antibody escape or receptor affinity selections
assays = antibody_escape_config["assays"]
assert "func_effects" not in assays, "cannot have assay called 'func_effects'"
assay_selections = {}
avg_assay_config = {}
for assay, assay_d in assays.items():
    assay_selections[assay] = antibody_escape_config[assay_d["selections"]]
    if "averages" in assay_d:
        avg_assay_config[assay] = antibody_escape_config[assay_d["averages"]]

#  make sure all samples defined
for sel_name, sel_d in itertools.chain.from_iterable(
    assay_selections[assay].items() for assay in assays
):
    for s in [sel_d["no_antibody_sample"], *list(sel_d["antibody_samples"])]:
        if s not in sample_to_library:
            raise ValueError(f"sample {s} for {selection_name} not in barcode_runs")

# Names and values of files to add to docs
assay_docs = {
    assay: {
        "Final summary plots": collections.defaultdict(dict),
        "Analysis notebooks": collections.defaultdict(dict),
        "Data files": collections.defaultdict(dict),
    }
    for assay in assay_selections
    if len(assay_selections[assay])
}


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
        nb=os.path.join(config["pipeline_path"], "notebooks/fit_escape.ipynb"),
    output:
        prob_escape_mean="results/{assay}/by_selection/{selection}_prob_escape_mean.csv",
        pickle="results/{assay}/by_selection/{selection}_polyclonal_model.pickle",
        nb="results/notebooks/fit_escape_{assay}_{selection}.ipynb",
    params:
        params_yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "params": assay_selections[wc.assay][wc.selection],
                "neut_standard_frac_csvs": list(input.neut_standard_fracs),
                "prob_escape_csvs": list(input.prob_escapes),
                "assay_config": assays[wc.assay],
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
            -p assay {wildcards.assay} \
            -p selection {wildcards.selection} \
            -y "{params.params_yaml}" \
            &> {log}
        """


for assay, sels in assay_selections.items():
    for sel in sels:
        for sample in sels[sel]["antibody_samples"]:
            assay_docs[assay]["Data files"][
                "Probability (fraction) escape for each variant in each "
                + assay.replace("_", " ")
                + " selection (CSVs)"
            ][f"{sel} {sample}"] = rules.prob_escape.output.prob_escape.format(
                assay=assay, selection=sel, sample=sample
            )
            assay_docs[assay]["Analysis notebooks"][
                "Fits of polyclonal models to individual "
                + assay.replace("_", " ")
                + " selections"
            ][sel] = rules.fit_escape.output.nb.format(assay=assay, selection=sel)


rule avg_escape:
    """Average antibody escape or receptor affinity across several selections."""
    input:
        plot_hide_stats=lambda wc: [
            d["csv"]
            for d in avg_assay_config[wc.assay][wc.antibody][
                "plot_hide_stats"
            ].values()
        ],
        prob_escape_means=lambda wc: [
            rules.fit_escape.output.prob_escape_mean.format(
                assay=wc.assay, selection=sel
            )
            for sel in avg_assay_config[wc.assay][wc.antibody]["selections"]
        ],
        pickles=lambda wc: [
            rules.fit_escape.output.pickle.format(assay=wc.assay, selection=sel)
            for sel in avg_assay_config[wc.assay][wc.antibody]["selections"]
        ],
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/avg_escape.ipynb"),
    output:
        pickle="results/{assay}/averages/{antibody}_polyclonal_model.pickle",
        effect_csv="results/{assay}/averages/{antibody}_mut_effect.csv",
        icXX_csv="results/{assay}/averages/{antibody}_mut_icXX.csv",
        effect_html="results/{assay}/averages/{antibody}_mut_effect_nolegend.html",
        icXX_html="results/{assay}/averages/{antibody}_mut_icXX_nolegend.html",
        nb="results/notebooks/avg_escape_{assay}_{antibody}.ipynb",
    params:
        params_yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "params": avg_assay_config[wc.assay][wc.antibody],
                "prob_escape_mean_csvs": list(input.prob_escape_means),
                "pickles": list(input.pickles),
                "assay_config": assays[wc.assay],
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/avg_escape_{assay}_{antibody}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p assay {wildcards.assay} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p avg_pickle_file {output.pickle} \
            -p effect_csv {output.effect_csv} \
            -p icXX_csv {output.icXX_csv} \
            -p effect_html {output.effect_html} \
            -p icXX_html {output.icXX_html} \
            -y '{params.params_yaml}' \
            &> {log}
        """


for assay in avg_assay_config:
    assay_str = assay.replace("_", " ")
    for heading, sec, fname, is_icXX in [
        (
            f"Average selections for {assay_str}",
            "Analysis notebooks",
            rules.avg_escape.output.nb,
            False,
        ),
        (f"{assay_str} CSVs", "Data files", rules.avg_escape.output.effect_csv, False),
        (
            f"{assay_str} ICXX CSVs",
            "Data files",
            rules.avg_escape.output.icXX_csv,
            True,
        ),
        (
            f"{assay_str} mutation effect plots",
            "Final summary plots",
            "results/{assay}/averages/{antibody}_mut_effect.html",
            False,
        ),
        (
            f"{assay_str} mutation ICXX plots",
            "Final summary plots",
            "results/{assay}/averages/{antibody}_mut_icXX.html",
            True,
        ),
    ]:
        for antibody, antibody_config in avg_assay_config[assay].items():
            if (not is_icXX) or (
                ("show_icXX_in_docs" in antibody_config)
                and antibody_config["show_icXX_in_docs"]
            ):
                assay_docs[assay][sec][heading][antibody] = fname.format(
                    assay=assay, antibody=antibody
                )

for assay, assay_doc in assay_docs.items():
    if assay_doc:
        docs[assays[assay]["title"]] = assay_doc
