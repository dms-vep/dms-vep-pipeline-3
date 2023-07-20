"""Implements ``snakemake`` rule to translate gene sequence."""


import sys

import Bio.SeqIO


sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

gene = Bio.SeqIO.read(snakemake.input.gene, "fasta").seq

# translate, making sure it gives valid protein with no stop except possibly at end
if len(gene) % 3 != 0:
    raise ValueError(f"{len(gene)} is not divisible by 3")
protseq = gene.translate()
assert len(protseq) == len(gene) // 3
if "*" in protseq[:-1]:
    raise ValueError(f"Premature stop codon in protein:\n{protseq}")

with open(snakemake.output.prot, "w") as f:
    f.write(f">gene\n{str(protseq)}\n")
