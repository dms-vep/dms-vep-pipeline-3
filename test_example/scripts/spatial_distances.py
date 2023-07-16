"""Implements ``snakemake`` rule `spatial_distances`"""


import sys
import urllib

import polyclonal.pdb_utils


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

print(f"Geting {snakemake.output.pdb=} from {snakemake.params.url=}")
urllib.request.urlretrieve(snakemake.params.url, snakemake.output.pdb)

print(f"Calculating distances using {snakemake.params.target_chains=}")
spatial_distances = polyclonal.pdb_utils.inter_residue_distances(
    snakemake.output.pdb,
    target_chains=snakemake.params.target_chains,
    target_atom=None,
)

print(f"Writing distances to {snakemake.output.csv=}")
spatial_distances.to_csv(snakemake.output.csv, index=False, float_format="%.5g")
