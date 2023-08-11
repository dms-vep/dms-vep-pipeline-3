"""``snakemake`` files with rules for calculating antibody/sera escape."""


# read the config for antibody escape
with open(config["antibody_escape_config"]) as f:
    antibody_escape_config = yaml.safe_load(f)

# get configuration for antibody escape and make sure all samples defined
antibody_selections = antibody_escape_config["antibody_selections"]
for selection_name, selection_d in antibody_selections.items():
    for s in [
        selection_d["no_antibody_sample"],
        *list(selection_d["antibody_samples"]),
    ]:
        if s not in sample_to_library:
            raise ValueError(f"sample {s} for {selection_name} not in barcode_runs")

# Names and values of files to add to docs
antibody_escape_docs = collections.defaultdict(dict)


rule prob_escape_antibody:
    """Compute probability (fraction) antibody escape for each variant."""
    input:
        no_antibody_sample=lambda wc: os.path.join(
            "results/barcode_counts",
            f"{antibody_selections[wc.selection]['no_antibody_sample']}_counts.csv",
        ),
        antibody_sample="results/barcode_counts/{sample}_counts.csv",
        codon_variants=config["codon_variants"],
        gene_sequence_codon=config["gene_sequence_codon"],
        site_numbering_map=config["site_numbering_map"],
    output:
        **{
            metric: f"results/antibody_escape/by_selection/{{selection}}/{{sample}}_{metric}.csv"
            for metric in ["prob_escape", "neut_standard_fracs"]
        },
    params:
        neut_standard=lambda wc: antibody_selections[wc.selection]["neut_standard_name"],
        # script checks that dates and libraries match for all samples
        dates=lambda wc: {
            "no_antibody_sample": sample_to_date[
                antibody_selections[wc.selection]["no_antibody_sample"]
            ],
            "antibody_sample": sample_to_date[wc.sample],
        },
        libraries=lambda wc: {
            "no_antibody_sample": sample_to_library[
                antibody_selections[wc.selection]["no_antibody_sample"]
            ],
            "antibody_sample": sample_to_library[wc.sample],
        },
    conda:
        "environment.yml"
    log:
        "results/logs/prob_escape_{selection}/{sample}.txt",
    script:
        "scripts/prob_escape.py"


for sel in antibody_selections:
    for sample in antibody_selections[sel]["antibody_samples"]:
        antibody_escape_docs[
            "Probability (fraction) escape for each variant (CSVs) in each selection"
        ][f"{sel} {sample}"] = rules.prob_escape_antibody.output.prob_escape.format(
            selection=sel, sample=sample
        )


rule fit_antibody_escape:
    """Fit a ``polyclonal`` model to an antibody selection."""
    input:
        spatial_distances=lambda wc: (
            [
                antibody_selections[wc.selection]["polyclonal_params"][
                    "spatial_distances"
                ]
            ]
            if antibody_selections[wc.selection]["polyclonal_params"][
                "spatial_distances"
            ]
            else []
        ),
        plot_hide_stats=lambda wc: [
            d["csv"]
            for d in antibody_selections[wc.selection]["plot_hide_stats"].values()
        ],
        prob_escapes=lambda wc: [
            rules.prob_escape_antibody.output.prob_escape.format(
                selection=wc.selection, sample=sample
            )
            for sample in antibody_selections[wc.selection]["antibody_samples"]
        ],
        neut_standard_fracs=lambda wc: [
            rules.prob_escape_antibody.output.neut_standard_fracs.format(
                selection=wc.selection, sample=sample
            )
            for sample in antibody_selections[wc.selection]["antibody_samples"]
        ],
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/fit_antibody_escape.ipynb"),
    output:
        prob_escape_mean="results/antibody_escape/by_selection/{selection}_prob_escape_mean.csv",
        pickle="results/antibody_escape/by_selection/{selection}_polyclonal_model.pickle",
        nb="results/notebooks/fit_antibody_escape_{selection}.ipynb",
    params:
        params_yaml=lambda wc: yaml.dump({"params": antibody_selections[wc.selection]}),
    conda:
        "environment.yml"
    log:
        "results/logs/fit_antibody_escape_{selection}.txt",
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


antibody_escape_docs["Fits of polyclonal models to antibody escape selections"] = {
    s: f"results/notebooks/fit_antibody_escape_{s}.ipynb" for s in antibody_selections
}

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
            rules.fit_antibody_escape.output.prob_escape_mean.format(selection=sel)
            for sel in avg_antibody_escape_config[wc.antibody]["selections"]
        ],
        pickles=lambda wc: [
            rules.fit_antibody_escape.output.pickle.format(selection=sel)
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
        antibody_escape_docs[heading][antibody] = fname.format(antibody=antibody)

docs["Antibody/serum escape"] = antibody_escape_docs
