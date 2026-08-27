"""Implements ``snakemake`` rule to translate gene sequence."""

import sys

import Bio.SeqIO

sys.stderr = sys.stdout = log = open(snakemake.log[0], "w")

gene = Bio.SeqIO.read(snakemake.input.gene, "fasta").seq

# codon positions (1-indexed) allowed to be a stop codon in addition to the
# last codon, e.g. for a gene spanning a natural readthrough stop codon
allowed_internal_stop_codons = set(getattr(snakemake.params, "allowed_internal_stop_codons", []))

# translate, making sure it gives valid protein with no unexpected stop codons
if len(gene) % 3 != 0:
    raise ValueError(f"{len(gene)} is not divisible by 3")
protseq = gene.translate()
assert len(protseq) == len(gene) // 3
unexpected_stops = sorted(
    i + 1
    for i, aa in enumerate(protseq[:-1])
    if aa == "*" and (i + 1) not in allowed_internal_stop_codons
)
if unexpected_stops:
    raise ValueError(
        f"Unexpected premature stop codon(s) at codon position(s) {unexpected_stops} "
        f"in protein:\n{protseq}"
    )

with open(snakemake.output.prot, "w") as f:
    f.write(f">gene\n{str(protseq)}\n")
