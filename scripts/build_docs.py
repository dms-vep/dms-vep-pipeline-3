"""Implements ``snakemake`` rule to translate gene sequence."""


import sys

import markdown


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

repo_url = snakemake.params.github_repo_url

md_text = [
    f"# {snakemake.params.description}",
    f"Analysis by {snakemake.params.authors} ({snakemake.params.year})",
    "",
    f"See [{repo_url}]({repo_url}) for full code.",
]

heading_depth = 1  # keys this deep are headings
init_heading = "##"  # first heading level is this in markdown


def process_docs(d, depth):
    """Recursive function to process ``docs_links`` nested dict."""
    depth += 1
    depth_diff = depth - heading_depth
    for key, val in d.items():
        if isinstance(val, str):
            entry = f" [{key}]({val})"
        elif isinstance(val, dict):
            entry = f" {key}"
        else:
            raise ValueError(f"{key=} has invalid value type {type(val)}\n{val}")
        if depth_diff <= 0:
            md_text.append(init_heading + "#" * (-depth_diff) + entry)
        else:
            md_text.append(" " * depth_diff + f"- {entry}")
        if isinstance(val, dict):
            process_docs(val, depth)


process_docs(snakemake.params.docs_links, 0)

md_text = "\n".join(md_text)

print(f"Rendering the following markdown text:\n\n{md_text}")

html = markdown.markdown(md_text)

with open(snakemake.output.html, "w") as f:
    f.write(html)
