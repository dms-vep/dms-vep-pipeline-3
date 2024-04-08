"""``snakemake`` rules for summarizing results across assays."""

# read the config for the summaries
with open(config["summaries_config"]) as f:
    summaries_config = yaml.safe_load(f)


rule summary:
    """Summary across different assays."""
    input:
        unpack(
            lambda wc: {
                f"antibody_escape {antibody}": rules.avg_escape.output.effect_csv.format(
                    assay="antibody_escape",
                    antibody=antibody,
                )
                for antibody_set_d in summaries_config[wc.summary][
                    "antibody_escape"
                ].values()
                for antibody in antibody_set_d["antibody_list"]
            }
            | {
                f"func_effects {condition_d['condition']}": os.path.join(
                    "results/func_effects/averages",
                        f"{condition_d['condition']}_{condition_d['effect_type']}.csv",
                    )
                    for condition_d in summaries_config[wc.summary][
                    "func_effects"
                ].values()
            }
            | {
                f"{assay} {condition_d['condition']}": rules.avg_escape.output.effect_csv.format(
                    assay=assay,
                    antibody=condition_d["condition"],
                )
                for (assay, assay_d) in summaries_config[wc.summary][
                    "other_assays"
                ].items()
                for condition_d in assay_d.values()
            }
        ),
        site_numbering_map_csv=config["site_numbering_map"],
        nb=os.path.join(config["pipeline_path"], "notebooks/summary.ipynb"),
    output:
        chart_overlaid="results/summaries/{summary}_overlaid_nolegend.html",
        chart_faceted="results/summaries/{summary}_faceted_nolegend.html",
        csv="results/summaries/{summary}.csv",
        per_antibody_escape_csv="results/summaries/{summary}_per_antibody_escape.csv",
        nb="results/notebooks/summary_{summary}.ipynb",
    params:
        yaml=lambda wc, input: yaml.round_trip_dump(
            {
                "config": {
                    key: val
                    for (key, val) in summaries_config[wc.summary].items()
                    if key not in ["title", "legend"]
                },
                "input_csvs": dict(input),
            }
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/summary_{summary}.txt",
    shell:
        """
        papermill {input.nb} {output.nb} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p chart_faceted {output.chart_faceted} \
            -p chart_overlaid {output.chart_overlaid} \
            -p output_csv_file {output.csv} \
            -p per_antibody_escape_csv {output.per_antibody_escape_csv} \
            -y "{params.yaml}" \
            &> {log}
        """


# Add files to docs
docs["Integrated summary plots"] = {
    summary.replace("_", " "): {
        **(
            {
                "Summary plot (escape overlaid)": f"results/summaries/{summary}_overlaid.html",
                "Summary plot (escape faceted)": f"results/summaries/{summary}_faceted.html",
            }
            if summaries_config[summary]["antibody_escape"]
            else {"Summary plot": f"results/summaries/{summary}_overlaid.html"}
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
