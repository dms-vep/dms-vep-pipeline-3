"""Implements ``snakemake`` rule to translate gene sequence."""


import sys

import Bio.SeqIO


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

gene = Bio.SeqIO.read(snakemake.input.gene, "fasta").seq

# first make sure gene can translate as CDS
_ = gene.translate(cds=True)
# now translate without cds flag to keep any stop codons
protseq = gene.translate()

with open(snakemake.output.prot, "w") as f:
    f.write(f">gene\n{str(protseq)}\n")
