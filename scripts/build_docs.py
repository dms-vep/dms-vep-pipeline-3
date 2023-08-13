"""Implements ``snakemake`` rule to translate gene sequence."""


import os
import sys

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
min_collapse_length = 1  # only collapse if at least this many entries


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
            if depth_diff > 0 and len(val) >= min_collapse_length:
                if depth_diff > 1:
                    # in order to handle larger list depths, the code to find
                    # `end_index` below would need to find the first matching </ul>
                    # rather than just the first one
                    raise ValueError("currently cannot handle list depth > 1")
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
        list_start = f"<li>{list_heading}<ul>"
        if html.count(list_start) != 1:
            raise ValueError(f"{html.count(list_start)} occurrences of {list_start}")
        # if list depth gets greater than one, need to start finding matching end
        end_index = (
            html[html.index(list_start) + len(list_start) :].index("</ul>")
            + html.index(list_start)
            + len(list_start)
        )

        to_replace = html[html.index(list_start) : end_index + len("</ul>")]
        assert html.count(to_replace) == 1
        assert to_replace.startswith("<li>")
        assert to_replace.endswith("</ul>"), to_replace

        replace_with = (
            "<li><details><summary>"
            + list_heading
            + " (click triangle to expand/collapse)</summary><ul>\n"
            + to_replace[len(list_start) :]
            + "</details>\n"
        )
        html = html.replace(to_replace, replace_with)

if os.path.dirname(snakemake.output.html):
    os.makedirs(os.path.dirname(snakemake.output.html), exist_ok=True)
with open(snakemake.output.html, "w") as f:
    f.write(html)
