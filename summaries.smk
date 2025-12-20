"""``snakemake`` rules for summarizing results across assays."""

# read the config for the summaries
with open(config["summaries_config"]) as f:
    summaries_config = {
        key: val
        for (key, val) in yaml.YAML(typ="safe", pure=True).load(f).items()
        if not key.endswith("_default")
    }


rule summary:
    """Summary across different assays."""
    input:
        **(
            {"mutation_annotations": config["mutation_annotations"]}
            if "mutation_annotations" in config
            else {}
        ),
        input_csvs=lambda wc: (
            [
                csv
                for key, val in summaries_config[wc.summary]["antibody_escape"].items()
                for csv in val["antibody_list"].values()
            ]
            + [
                val["csv"]
                for key, val in summaries_config[wc.summary][
                    "other_phenotypes"
                ].items()
            ]
        ),
        site_numbering_map=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/summary.ipynb"),
    output:
        chart_overlaid="results/summaries/{summary}_overlaid.html",
        chart_faceted="results/summaries/{summary}_faceted.html",
        csv="results/summaries/{summary}.csv",
        per_antibody_escape_csv="results/summaries/{summary}_per_antibody_escape.csv",
        nb="results/notebooks/summary_{summary}.ipynb",
    params:
        config_params_yaml=lambda wc: yaml_str(
            {
                "config": summaries_config[wc.summary],
            }
        ),
        input_params_yaml=lambda wc, input: yaml_str(
            {
                "input_config": {
                    key: list(val) if key == "input_csvs" else val
                    for key, val in dict(input).items()
                }
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/summary_{summary}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p chart_faceted {output.chart_faceted} \
            -p chart_overlaid {output.chart_overlaid} \
            -p output_csv_file {output.csv} \
            -p per_antibody_escape_csv {output.per_antibody_escape_csv} \
            -y "{params.config_params_yaml}" \
            -y "{params.input_params_yaml}" \
            &> {log}
        """


# Add files to docs
docs["Integrated summary plots"] = {
    summary.replace("_", " "): {
        **(
            {
                "Summary plot (escape overlaid)": rules.summary.output.chart_overlaid.format(
                    summary=summary
                ),
                "Summary plot (escape faceted)": rules.summary.output.chart_faceted.format(
                    summary=summary
                ),
            }
            if summaries_config[summary]["antibody_escape"]
            else {
                "Summary plot": rules.summary.output.chart_overlaid.format(
                    summary=summary
                )
            }
        ),
        "CSV summarizing results": rules.summary.output.csv.format(summary=summary),
        **(
            {
                "CSV summarizing per-antibody results": rules.summary.output.per_antibody_escape_csv.format(
                    summary=summary
                )
            }
            if summaries_config[summary]["antibody_escape"]
            else {}
        ),
    }
    for summary in summaries_config
}
