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
            metric: f"results/antibody_escape/{{selection}}/{{sample}}_{metric}.csv"
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
        "results/logs/prob_escape_{selection}_{sample}.txt",
    script:
        "scripts/prob_escape.py"


for sel in antibody_selections:
    for sample in antibody_selections[sel]["antibody_samples"]:
        antibody_escape_docs[
            "Probability (fraction) antibody escape for each variant (CSVs)"
        ][f"{sel} {sample}"] = f"results/antibody_escape/{sel}/{sample}_prob_escape.csv"


rule fit_antibody_escape:
    """Fit a ``polyclonal`` model to an antibody selection."""
    input:
        prob_escapes=lambda wc: [
            f"results/antibody_escape/{wc.selection}/{sample}_prob_escape.csv"
            for sample in antibody_selections[wc.selection]["antibody_samples"]
        ],
        neut_standard_fracs=lambda wc: [
            f"results/antibody_escape/{wc.selection}/{sample}_neut_standard_fracs.csv"
            for sample in antibody_selections[wc.selection]["antibody_samples"]
        ],
        nb=os.path.join(config["pipeline_path"], "notebooks/fit_antibody_escape.ipynb"),
    output:
        prob_escape_mean="results/antibody_escape/{selection}/prob_escape_mean.csv",
        nb="results/notebooks/fit_antibody_escape_{selection}.ipynb",
    params:
        params_yaml=lambda wc: yaml.dump({"params": antibody_selections[wc.selection]}),
    conda:
        "environment.yml"
    log:
        "results/fit_antibody_escape_{selection}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p prob_escape_mean_csv {output.prob_escape_mean} \
            -p selection {wildcards.selection} \
            -y "{params.params_yaml}" \
            &> {log}
        """


antibody_escape_docs["Fits of polyclonal models to antibody escape selections"] = {
    s: f"results/notebooks/fit_antibody_escape_{s}.ipynb" for s in antibody_selections
}


docs["Antibody/serum escape"] = antibody_escape_docs
