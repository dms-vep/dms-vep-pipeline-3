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
    """Convert a Jupyter notebook to HTML in its results directory."""
    input:
        nb="results/notebooks/{notebook}.ipynb",
    output:
        html="results/notebooks/{notebook}.html",
    conda:
        "environment.yml"
    log:
        "results/logs/notebook_html_{notebook}.txt",
    shell:
        "jupyter nbconvert --to html {input.nb} &> {log}"


rule build_docs:
    """Build the HTML documentation."""
    input:
        docs_processed_files.values(),
    output:
        docs_dir=directory(config["docs"]),
        html=os.path.join(config["docs"], "index.html"),
    params:
        github_repo_url=config["github_repo_url"],
        docs_links=docs_links,
        docs_processed_files=docs_processed_files,
        description=config["description"],
        year=config["year"],
        authors=config["authors"],
    conda:
        "environment.yml"
    log:
        "results/logs/build_docs.txt",
    script:
        "scripts/build_docs.py"


if "build_vitepress_homepage" in config and config["build_vitepress_homepage"]:

    rule build_vitepress_homepage:
        """Copy all files from the docs directory to the VitePress homepage directory"""
        input:
            html=os.path.join(config["docs"], "index.html"),
        output:
            html=os.path.join(config["homepage"], "appendix.html"),
        params:
            docs=lambda _, input: os.path.dirname(input.html),
            homepage=lambda _, output: os.path.dirname(output.html),
        log:
            "results/logs/build_vitepress_homepage.txt",
        conda:
            "environment.yml"
        shell:
            """
            # Copy contents of docs/ to homepage/public/
            cp -r {params.docs}/* {params.homepage}
            # Remove the index.html file
            rm -f {params.homepage}/index.html
            # Copy and rename the index.html file
            cp {input.html} {output.html}
            """

    other_target_files.append(rules.build_vitepress_homepage.output.html)
