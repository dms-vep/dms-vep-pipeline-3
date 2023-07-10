"""Implements ``snakemake`` rule to translate gene sequence."""


import sys

import bs4

import markdown
import markdown.extensions.toc


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

repo_url = snakemake.params.github_repo_url

md_text = [
    f"# {snakemake.params.description}",
    f"Analysis by {snakemake.params.authors} ({snakemake.params.year})",
    "",
    f"See [{repo_url}]({repo_url}) for full code.",
    "",
    # table of contents: https://python-markdown.github.io/extensions/toc/
    "[TOC]",
    "",
]

subheading_depth = 1  # keys this deep are subheadings
init_subheading = "##"  # first subheading level is this in markdown
collapse_list = []  # nested list headings to collapse


def process_docs(d, depth):
    """Recursive function to process ``docs_links`` nested dict."""
    depth += 1
    depth_diff = depth - subheading_depth
    for key, val in d.items():
        if isinstance(val, str):
            entry = f" [{key}]({val})"
        elif isinstance(val, dict):
            entry = f" {key}"
        else:
            raise ValueError(f"{key=} has invalid value type {type(val)}\n{val}")
        if depth_diff <= 0:
            md_text.append(
                init_subheading + "#" * (subheading_depth + depth_diff - 1) + entry
            )
        else:
            md_text.append("  " * depth_diff + f"-{entry}")
        if isinstance(val, dict):
            process_docs(val, depth)
            if depth_diff > 0:
                collapse_list.append(key)


process_docs(snakemake.params.docs_links, 0)

md_text = "\n".join(md_text)

print(f"Rendering the following markdown text:\n\n{md_text}\n\n")

html = markdown.markdown(
    md_text,
    extensions=[
        markdown.extensions.toc.TocExtension(
            title="Contents",
            toc_depth="2-4",
        ),
    ],
)

# edit HTML to make deeply nested list items collapsible:
# https://gist.github.com/dotiful/0bd3516f42c6ca68479e64ad2942ac90
collapse_nested_lists = True
if collapse_nested_lists:
    print(f"Collapsing the following nested lists:\n{collapse_list}")
    for list_heading in collapse_list:
        tags = [str(tag) for tag in bs4.BeautifulSoup(html, "html.parser")]
        to_collapse = [
            tag for tag in tags if tag.startswith(f"<ul>\n<li>{list_heading}<ul>")
        ]
        if len(to_collapse) != 1:
            raise ValueError(f"not 1 head for {list_heading}\n{to_collapse=}\n{tags=}")
        to_collapse = to_collapse[0]
        start = f"<ul>\n<li>{list_heading}<ul>\n"
        end = "\n</ul>"
        assert to_collapse.startswith(start)
        assert to_collapse.endswith(end)
        collapsed = (
            "<ul>\n<li><details><summary>"
            + list_heading
            + " (click triangle to expand/collapse)</summary><ul>\n"
            + to_collapse[len(start) : -len(end)]
            + "</details>\n</ul>"
        )
        assert html.count(to_collapse) == 1
        html = html.replace(to_collapse, collapsed)

with open(snakemake.output.html, "w") as f:
    f.write(html)
