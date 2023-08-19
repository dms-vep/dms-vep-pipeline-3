"""Functions for ``docs.smk``."""


def process_nested_docs_dict(d, github_blob_url):
    """Recursive function to process docs dict to get input files and links.

    Parameters
    ----------
    d : dict
        Nested dictionary containing headings and entries for documentation.
    github_blob_url : str
        Base path to files stored on GitHub for this repo.

    Returns
    -------
    (d_links, processed_files)
        `d_links` is a copy of `d` where final values have appropriate paths for docs.
        `processed_files` is a dict keyed by names of files as they are processed and
        copied to docs, with values original file name.

    """
    d_links = {}
    processed_files = {}
    for key, val in d.items():
        if isinstance(val, str):
            path, ext = os.path.splitext(val)
            if ext == ".gz":
                gz = True
                path, ext = os.path.splitext(path)
            else:
                gz = False
            base = os.path.basename(path)
            processed_f = None
            if ext == ".ipynb" and not gz:
                d_links[key] = f"notebooks/{base}.html"
                processed_f = os.path.join(config["docs"], d_links[key])
            elif ext == ".html" and not gz:
                d_links[key] = f"htmls/{base}.html"
                processed_f = os.path.join(config["docs"], d_links[key])
            elif ext in [".csv", ".fasta", ".fa"]:
                d_links[key] = os.path.join(github_blob_url, val)
            else:
                raise ValueError(
                    f"cannot handle file extension {ext=} and {gz=} as for {val=}"
                )
            if processed_f is not None:
                assert (
                    processed_f not in processed_files
                ), f"duplicate {processed_f}\n{d}"
                processed_files[processed_f] = val
        elif isinstance(val, dict):
            d_links[key], pfiles = process_nested_docs_dict(val, github_blob_url)
            dup_files = set(processed_files).intersection(pfiles)
            if dup_files:
                raise f"duplicate processed file names {dup_files}"
            processed_files.update(pfiles)
        else:
            raise ValueError(f"value for {key=} is invalid type {type(val)}\n{val=}")

    return d_links, processed_files
