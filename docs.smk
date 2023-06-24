"""``snakemake`` files with rules for building HTML documentation.

These rules operation on the nested dict ``docs``, which specifies
the structure of the documentation to build.

"""


# functions used by these rulesj
include: "docs_funcs.smk"


docs_links, docs_processed_files = process_nested_docs_dict(
    docs,
    config["github_blob_url"],
)


rule notebook_html:
    """Convert a Jupyter notebook into a HTML for docs."""
    input:
        nb=lambda wc: docs_processed_files[
            os.path.join(config["docs"], f"notebooks/{wc.notebook}.html")
        ],
    output:
        html=os.path.join(config["docs"], "notebooks/{notebook}.html"),
    params:
        subdir=lambda _, output: os.path.dirname(output.html),
    conda:
        "environment.yml"
    log:
        "results/logs/notebook_html_{notebook}.txt",
    shell:
        "jupyter nbconvert --to html --output-dir {params.subdir} {input.nb} &> {log}"


rule cp_html:
    """Copy a HTML file for docs."""
    input:
        html=lambda wc: docs_processed_files[
            os.path.join(config["docs"], f"htmls/{wc.html}.html")
        ],
    output:
        html=os.path.join(config["docs"], "htmls/{html}.html"),
    conda:
        "environment.yml"
    log:
        "results/logs/cp_html_{html}.txt",
    shell:
        "cp {input.html} {output.html}"


rule build_docs:
    """Build the HTML documentation."""
    input:
        flatdict.FlatDict(docs).values(),
        docs_processed_files.keys(),
    output:
        html=os.path.join(config["docs"], "index.html"),
    params:
        github_repo_url=config["github_repo_url"],
        docs_links=docs_links,
        description=config["description"],
        year=config["year"],
        authors=config["authors"],
    conda:
        "environment.yml"
    log:
        "results/logs/build_docs.txt",
    script:
        "scripts/build_docs.py"
