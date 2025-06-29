"""``snakemake`` files with rules for calculating functional effects."""

# read the config for func effects
with open(config["func_effects_config"]) as f:
    func_effects_config = yaml.YAML(typ="safe", pure=True).load(f)

# get configuration for functional scores and make sure all samples defined
func_scores = func_effects_config["func_scores"]
for selection_name, selection_d in func_scores.items():
    for s in ["post_selection_sample", "pre_selection_sample"]:
        if selection_d[s] not in sample_to_library:
            raise ValueError(
                f"{s} for {selection_name} of {selection_d[s]} not in barcode_runs"
            )

# Names and values of files to add to docs
func_effects_docs = {
    "Final summary plots": collections.defaultdict(dict),
    "Analysis notebooks": collections.defaultdict(dict),
    "Data files": collections.defaultdict(dict),
}


rule func_scores:
    """Compute functional scores for variants."""
    input:
        unpack(
            lambda wc: {
                s: f"results/barcode_counts/{func_scores[wc.selection][s]}_counts.csv"
                for s in ["post_selection_sample", "pre_selection_sample"]
            }
        ),
        codon_variants=config["codon_variants"],
        gene_sequence_codon=config["gene_sequence_codon"],
        site_numbering_map=config["site_numbering_map"],
    output:
        func_scores="results/func_scores/{selection}_func_scores.csv",
        count_summary="results/func_scores/{selection}_count_summary.csv",
    params:
        func_score_params=lambda wc: func_scores[wc.selection]["func_score_params"],
        samples=lambda wc: {
            s: func_scores[wc.selection][s]
            for s in ["post_selection_sample", "pre_selection_sample"]
        },
        dates=lambda wc: {
            s: sample_to_date[func_scores[wc.selection][s]]
            for s in ["post_selection_sample", "pre_selection_sample"]
        },
        # script will throw error if pre_library and post_library differ
        libraries=lambda wc: {
            s: sample_to_library[func_scores[wc.selection][s]]
            for s in ["post_selection_sample", "pre_selection_sample"]
        },
    conda:
        "environment.yml"
    log:
        "results/logs/func_scores_{selection}.txt",
    script:
        "scripts/func_scores.py"


for s in func_scores:
    func_effects_docs["Data files"]["Per-variant functional score CSVs"][
        s
    ] = f"results/func_scores/{s}_func_scores.csv"


rule analyze_func_scores:
    """Analyze functional scores."""
    input:
        func_scores=expand(rules.func_scores.output.func_scores, selection=func_scores),
        count_summaries=expand(
            rules.func_scores.output.count_summary,
            selection=func_scores,
        ),
        nb=os.path.join(config["pipeline_path"], "notebooks/analyze_func_scores.ipynb"),
    output:
        nb="results/notebooks/analyze_func_scores.ipynb",
    params:
        func_scores=yaml_str({"selections": func_effects_config["func_scores"]}),
    conda:
        "environment.yml"
    log:
        "results/logs/analyze_func_scores.txt",
    shell:
        """
        papermill {input.nb} {output.nb} -y "{params.func_scores}" &> {log}
        """


func_effects_docs["Analysis notebooks"][
    "Analysis of functional scores"
] = rules.analyze_func_scores.output.nb


rule func_effects_global_epistasis:
    """Fit global epistasis model to func scores to get mutation functional effects."""
    input:
        func_scores="results/func_scores/{selection}_func_scores.csv",
        site_numbering_map=config["site_numbering_map"],
        nb=os.path.join(
            config["pipeline_path"],
            "notebooks/func_effects_global_epistasis.ipynb",
        ),
    output:
        func_effects="results/func_effects/by_selection/{selection}_func_effects.csv",
        nb="results/notebooks/func_effects_global_epistasis_{selection}.ipynb",
    params:
        global_epistasis_params_yaml=lambda wc: yaml_str(
            {
                "global_epistasis_params": func_scores[wc.selection][
                    "global_epistasis_params"
                ],
            }
        ),
    threads: 1
    conda:
        "environment.yml"
    log:
        "results/logs/func_effects_global_epistasis_{selection}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p selection {wildcards.selection} \
            -p func_scores {input.func_scores} \
            -p func_effects {output.func_effects} \
            -p site_numbering_map {input.site_numbering_map} \
            -p threads {threads} \
            -y "{params.global_epistasis_params_yaml}" \
            &> {log}
        """


for s in func_scores:
    func_effects_docs["Analysis notebooks"][
        "Notebooks with per-selection global epistasis fitting"
    ][s] = f"results/notebooks/func_effects_global_epistasis_{s}.ipynb"
    func_effects_docs["Data files"]["Per-selection mutation functional effect CSVs"][
        s
    ] = f"results/func_effects/by_selection/{s}_func_effects.csv"


rule avg_func_effects:
    """Average and plot the functional effects for a condition."""
    input:
        **(
            {"mutation_annotations_csv": config["mutation_annotations"]}
            if "mutation_annotations" in config
            else {}
        ),
        site_numbering_map_csv=config["site_numbering_map"],
        selections=lambda wc: [
            f"results/func_effects/by_selection/{s}_func_effects.csv"
            for s in func_effects_config["avg_func_effects"][wc.condition][
                "selections"
            ]
        ],
        nb=os.path.join(config["pipeline_path"], "notebooks/avg_func_effects.ipynb"),
    output:
        nb="results/notebooks/avg_func_effects_{condition}.ipynb",
        func_effects_csv="results/func_effects/averages/{condition}_func_effects.csv",
        func_effects_singlemut_csv="results/func_effects/averages/{condition}_func_effects_singlemut.csv",
        latent_effects_csv="results/func_effects/averages/{condition}_latent_effects.csv",
        functional_html="results/func_effects/averages/{condition}_func_effects.html",
        functional_singlemut_html="results/func_effects/averages/{condition}_func_effects_singlemut.html",
        latent_html="results/func_effects/averages/{condition}_latent_effects.html",
    params:
        params_yaml=lambda wc, input: yaml_str(
            {
                "params": func_effects_config["avg_func_effects"][wc.condition],
                "mutation_annotations_csv": (
                    input.mutation_annotations_csv
                    if input.get("mutation_annotations_csv")
                    else None
                ),
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/avg_func_effects_{condition}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p func_effects_csv {output.func_effects_csv} \
            -p func_effects_singlemut_csv {output.func_effects_singlemut_csv} \
            -p latent_effects_csv {output.latent_effects_csv} \
            -p functional_html {output.functional_html} \
            -p functional_singlemut_html {output.functional_singlemut_html} \
            -p latent_html {output.latent_html} \
            -y '{params.params_yaml}' \
            &> {log}
        """


func_effects_docs["Analysis notebooks"][
    "Notebooks averaging mutation functional effects"
] = {
    c: f"results/notebooks/avg_func_effects_{c}.ipynb"
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Data files"]["Average mutation functional effects (CSV files)"] = {
    c: f"results/func_effects/averages/{c}_func_effects.csv"
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Data files"][
    "Average mutation functional effects from single-mutant variants only (CSV files)"
] = {
    c: f"results/func_effects/averages/{c}_func_effects_singlemut.csv"
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Data files"][
    "Average mutation latent-phenotype effects (CSV files)"
] = {
    c: f"results/func_effects/averages/{c}_func_effects.csv"
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Final summary plots"][
    "Interactive plots of average mutation functional effects"
] = {
    c: rules.avg_func_effects.output.functional_html.format(condition=c)
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Final summary plots"][
    "Interactive plots of average mutation functional effects from single-mutant variants only"
] = {
    c: rules.avg_func_effects.output.functional_singlemut_html.format(condition=c)
    for c in func_effects_config["avg_func_effects"]
}

func_effects_docs["Final summary plots"][
    "Interactive plots of average mutation latent-phenotype effects"
] = {
    c: rules.avg_func_effects.output.latent_html.format(condition=c)
    for c in func_effects_config["avg_func_effects"]
}

# are we doing func_effect_diffs?
if ("func_effect_diffs" in func_effects_config) and (
    func_effects_config["func_effect_diffs"] is not None
):
    func_effect_diffs = func_effects_config["func_effect_diffs"]
else:
    func_effect_diffs = {}


rule func_effect_diffs:
    """Simple differences in functional effects between two conditions."""
    input:
        lambda wc: [
            rules.func_effects_global_epistasis.output.func_effects.format(
                selection=sel
            )
            for c in ["condition_1", "condition_2"]
            for sel in func_effect_diffs[wc.comparison][c]["selections"]
        ],
        **(
            {"mutation_annotations_csv": config["mutation_annotations"]}
            if "mutation_annotations" in config
            else {}
        ),
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/func_effect_diffs.ipynb"),
    output:
        diffs="results/func_effect_diffs/{comparison}_diffs.csv",
        chart="results/func_effect_diffs/{comparison}_diffs.html",
        corr_chart="results/func_effect_diffs/{comparison}_diffs_corr.html",
        nb="results/notebooks/func_effect_diffs_{comparison}.ipynb",
    params:
        params_yaml=lambda wc, input: yaml_str(
            {
                "params": {"singlemut": False} | func_effect_diffs[wc.comparison],
                "mutation_annotations_csv": (
                    input.mutation_annotations_csv
                    if input.get("mutation_annotations_csv")
                    else None
                ),
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/func_effect_diffs_{comparison}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p diffs_csv {output.diffs} \
            -p chart_html {output.chart} \
            -p corr_chart_html {output.corr_chart} \
            -y '{params.params_yaml}' \
            &> {log}
        """


if func_effect_diffs:
    func_effects_docs["Analysis notebooks"][
        "Notebooks taking differences in functional effects between conditions"
    ] = {
        c: rules.func_effect_diffs.output.nb.format(comparison=c)
        for c in func_effect_diffs
    }
    func_effects_docs["Data files"]["Differences in functional effects (CSV files)"] = {
        c: rules.func_effect_diffs.output.diffs.format(comparison=c)
        for c in func_effect_diffs
    }
    func_effects_docs["Final summary plots"][
        "Interactive plots of differences in functional effects between conditions"
    ] = {
        c: rules.func_effect_diffs.output.chart.format(comparison=c)
        for c in func_effect_diffs
    } | {
        f"{c} correlation chart": rules.func_effect_diffs.output.corr_chart.format(
            comparison=c
        )
        for c in func_effect_diffs
    }


# are we doing func_effect_shifts comparisons?
if ("func_effect_shifts" in func_effects_config) and (
    func_effects_config["func_effect_shifts"] is not None
):
    func_effect_shifts = func_effects_config["func_effect_shifts"]
    if ("avg_func_effect_shifts" in func_effects_config) and (
        func_effects_config["avg_func_effect_shifts"] is not None
    ):
        avg_func_effect_shifts = func_effects_config["avg_func_effect_shifts"]
    else:
        avg_func_effect_shifts = {}
else:
    func_effect_shifts = {}
    avg_func_effect_shifts = {}


rule func_effect_shifts:
    """``multidms`` comparison of conditions to get shifts in functional effects."""
    input:
        lambda wc: [
            rules.func_scores.output.func_scores.format(selection=sel)
            for sel in func_effect_shifts[wc.comparison]["conditions"].values()
        ],
        site_numbering_map=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/func_effect_shifts.ipynb"),
    output:
        shifts="results/func_effect_shifts/by_comparison/{comparison}_shifts.csv",
        nb="results/notebooks/func_effect_shifts_{comparison}.ipynb",
    params:
        params_yaml=lambda wc: yaml_str({"params": func_effect_shifts[wc.comparison]}),
    threads: 1
    conda:
        "environment.yml"
    log:
        "results/logs/func_effect_shifts_{comparison}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -y "{params.params_yaml}" \
            -p site_numbering_map {input.site_numbering_map} \
            -p shifts_csv {output.shifts} \
            -p threads {threads} \
            &> {log}
        """


if func_effect_shifts:
    func_effects_docs["Analysis notebooks"][
        "Notebooks fitting shifts in functional effects"
    ] = {
        c: rules.func_effect_shifts.output.nb.format(comparison=c)
        for c in func_effect_shifts
    }
    func_effects_docs["Data files"]["Per-condition functional effect shifts CSVs"] = {
        c: rules.func_effect_shifts.output.shifts.format(comparison=c)
        for c in func_effect_shifts
    }


rule avg_func_effect_shifts:
    """Average and plot the functional effects shifts for a comparison."""
    input:
        lambda wc: [
            rules.func_effect_shifts.output.shifts.format(comparison=c)
            for c in avg_func_effect_shifts[wc.comparison]["comparisons"]
        ],
        **(
            {"mutation_annotations_csv": config["mutation_annotations"]}
            if "mutation_annotations" in config
            else {}
        ),
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(
            config["pipeline_path"],
            "notebooks/avg_func_effect_shifts.ipynb",
        ),
    output:
        shifts_csv="results/func_effect_shifts/averages/{comparison}_shifts.csv",
        shifts_html="results/func_effect_shifts/averages/{comparison}_shifts.html",
        nb="results/notebooks/avg_func_effect_shifts_{comparison}.ipynb",
    params:
        params_yaml=lambda wc, input: yaml_str(
            {
                "params": avg_func_effect_shifts[wc.comparison],
                "mutation_annotations_csv": (
                    input.mutation_annotations_csv
                    if input.get("mutation_annotations_csv")
                    else None
                ),
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/avg_func_effect_shifts_{comparison}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p shifts_csv {output.shifts_csv} \
            -p shifts_html {output.shifts_html} \
            -y '{params.params_yaml}' \
            &> {log}
        """


if avg_func_effect_shifts:
    func_effects_docs["Analysis notebooks"][
        "Notebooks averaging shifts in functional effects"
    ] = {
        c: rules.avg_func_effect_shifts.output.nb.format(comparison=c)
        for c in avg_func_effect_shifts
    }
    func_effects_docs["Data files"][
        "Average shifts in functional effects (CSV files)"
    ] = {
        c: rules.avg_func_effect_shifts.output.shifts_csv.format(comparison=c)
        for c in avg_func_effect_shifts
    }
    func_effects_docs["Final summary plots"][
        "Interactive plots of average shifts in functional effects"
    ] = {
        c: rules.avg_func_effect_shifts.output.shifts_html.format(comparison=c)
        for c in avg_func_effect_shifts
    }


docs["Functional effects of mutations"] = func_effects_docs
