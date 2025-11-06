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
        html="results/docs/index.html",
    params:
        github_repo_url=config["github_repo_url"],
        docs_links=docs_links,
        docs_processed_files=docs_processed_files,
        description=config["description"],
        year=config["year"],
        authors=config["authors"],
        docs_dir=lambda _, output: os.path.dirname(output.html),
    conda:
        "environment.yml"
    log:
        "results/logs/build_docs.txt",
    script:
        "scripts/build_docs.py"


if build_vitepress_homepage:

    rule build_vitepress_homepage:
        """Copy all files from the docs directory to the VitePress homepage directory"""
        input:
            html=rules.build_docs.output.html,
            homepage="homepage",
        output:
            html="results/homepage/public/appendix.html",
        params:
            docs=lambda _, input: os.path.dirname(input.html),
            homepage_public=lambda _, output: os.path.dirname(output.html),
            homepage=lambda _w, output: os.path.dirname(os.path.dirname(output.html)),
        log:
            "results/logs/build_vitepress_homepage.txt",
        conda:
            "environment.yml"
        shell:
            """
            rm -rf {params.homepage}
            cp -r {input.homepage} {params.homepage}
            cp -r {params.docs}/* {params.homepage_public}
            rm -f {params.homepage_public}/index.html
            cp {input.html} {output.html}
            """

    other_target_files.append(rules.build_vitepress_homepage.output.html)


rule build_publish_docs:
    """Get the actual docs to publish."""
    input:
        html=(
            rules.build_vitepress_homepage.output.html
            if build_vitepress_homepage
            else rules.build_docs.output.html
        ),
    output:
        publish_docs=directory("results/publish_docs"),
    params:
        input_dir=lambda _, input: (
            os.path.dirname(os.path.dirname(input.html))
            if build_vitepress_homepage
            else os.path.dirname(input.html)
        ),
    log:
        "results/logs/build_publish_docs.log",
    conda:
        "environment.yml"
    shell:
        """
        rm -rf {output.publish_docs} &> {log}
        cp -r {params.input_dir} {output.publish_docs} &> {log}
        """
