"""``snakemake`` rules for summarizing results across assays."""

# read the config for the summaries
with open(config["summaries_config"]) as f:
    summaries_config = yaml.safe_load(f)

summary_config = summaries_config["summary"]


rule summary:
    """Summary across different assays."""
    input:
        **{
            f"antibody_escape {antibody}": rules.avg_escape.output.effect_csv.format(
                assay="antibody_escape", antibody=antibody,
            )
            for antibody_set_d in summary_config["antibody_escape"].values()
            for antibody in antibody_set_d["antibody_list"]
        },
        **{
            f"func_effects {condition_d['condition']}": os.path.join(
                "results/func_effects/averages",
                condition_d["condition"] + "_" + condition_d["effect_type"] + ".csv",
            )
            for condition_d in summary_config["func_effects"].values()
        },
        **{
            f"{assay} {condition_d['condition']}": rules.avg_escape.output.effect_csv.format(
                assay=assay, antibody=condition_d["condition"],
            )
            for (assay, assay_d) in summary_config["other_assays"].items()
            for condition_d in assay_d.values()
        },
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/summary.ipynb"),
    output:
        chart_overlaid="results/summaries/summary_overlaid_nolegend.html",
        chart_faceted="results/summaries/summary_faceted_nolegend.html",
        csv="results/summaries/summary.csv",
        nb="results/notebooks/summary.ipynb",
    params:
        yaml=lambda _, input: yaml.round_trip_dump(
            {"config": summary_config, "input_csvs": dict(input)}
        )
    conda:
        "environment.yml"
    log:
        "results/logs/summary.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p chart_faceted {output.chart_faceted} \
            -p chart_overlaid {output.chart_overlaid} \
            -p output_csv_file {output.csv} \
            -y "{params.yaml}" \
            &> {log}
        """


# Add files to docs
docs["Summary of results across assays"] = {
    "Final summary plots": ( 
        {
            "Summary of assays (escape overlaid)": rules.summary.output.chart_overlaid,
            "Summary of assays (escape faceted)": rules.summary.output.chart_faceted,
        }
        if summary_config["antibody_escape"]
        else {"Summary of assays": rules.summary.output.chart_overlaid}
    ),
    "Data files": {
        "CSV summarizing results": rules.summary.output.csv,
    },
}