"""Custom rules used in the ``snakemake`` pipeline.

This file is included by the pipeline ``Snakefile``.

"""


rule spatial_distances:
    """Get spatial distances from PDB."""
    output:
        csv="results/spatial_distances/7tov.csv",
        pdb="results/spatial_distances/7tov.pdb",
    params:
        url=lambda _, output: os.path.join(
            "https://files.rcsb.org/download",
            os.path.basename(output.pdb),
        ),
        target_chains=["A", "B", "C"],
    log:
        log="results/logs/spatial_distances.txt",
    conda:
        os.path.join(config["pipeline_path"], "environment.yml")
    script:
        "scripts/spatial_distances.py"
