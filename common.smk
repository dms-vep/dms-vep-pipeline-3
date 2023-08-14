"""Common rules used in multiple ``snakemake`` ``*.smk`` files in pipeline."""


include: "common_funcs.smk"


rule format_altair_html:
    """Format ``altair`` charts by adding title, legend, etc."""
    input:
        html="results/{chart}_nolegend.html",
        pyscript=os.path.join(config["pipeline_path"], "scripts/format_altair_html.py"),
    output:
        html="results/{chart}.html",
        legend=temp("results/{chart}.md"),
    params:
        chart_params=format_altair_html_chart_params,
        suffix=(
            f"Analysis by {config['authors']} ({config['year']}).\n\n See "
            f"[{config['github_repo_url']}]({config['github_repo_url']}) for code/data."
        ),
    conda:
        "environment.yml"
    log:
        "results/logs/{chart}.txt",
    shell:
        """
        echo "## {params.chart_params[title]}\n" > {output.legend}
        echo "{params.chart_params[legend]}\n\n" >> {output.legend}
        echo "{params.suffix}" >> {output.legend}
        python {input.pyscript} \
            --chart {input.html} \
            --markdown {output.legend} \
            --title "{params.chart_params[title]}" \
            --output {output.html}
        """
