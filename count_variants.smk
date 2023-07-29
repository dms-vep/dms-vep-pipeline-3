"""``snakemake`` files with rules for counting variants from barcode sequencing."""


# Names and values of files to add to docs
count_variants_docs = collections.defaultdict(dict)


rule count_barcodes:
    """Count barcodes for each sample."""
    input:
        fastq_R1=lambda wc: barcode_runs.set_index("sample").at[wc.sample, "fastq_R1"],
        variants=config["codon_variants"],
    output:
        counts="results/barcode_counts/{sample}_counts.csv",
        invalid="results/barcode_counts/{sample}_invalid.csv",
        fates="results/barcode_counts/{sample}_fates.csv",
    params:
        parser_params=config["illumina_barcode_parser_params"],
        library=lambda wc: sample_to_library[wc.sample],
    conda:
        "environment.yml"
    log:
        "results/logs/count_barcodes_{sample}.txt",
    script:
        "scripts/count_barcodes.py"


for sample in barcode_runs["sample"]:
    count_variants_docs["Barcode counts"][sample] = os.path.join(
        "results/barcode_counts",
        f"{sample}_counts.csv",
    )


rule analyze_variant_counts:
    """Analysis of the counts of variants."""
    input:
        expand(rules.count_barcodes.output.counts, sample=barcode_runs["sample"]),
        expand(rules.count_barcodes.output.invalid, sample=barcode_runs["sample"]),
        expand(rules.count_barcodes.output.fates, sample=barcode_runs["sample"]),
        config["gene_sequence_codon"],
        config["codon_variants"],
        config["site_numbering_map"],
        config["barcode_runs"],
        nb=os.path.join(
            config["pipeline_path"],
            "notebooks/analyze_variant_counts.ipynb",
        ),
    output:
        nb="results/notebooks/analyze_variant_counts.ipynb",
    conda:
        "environment.yml"
    log:
        "results/logs/analyze_variant_counts.txt",
    shell:
        "papermill {input.nb} {output.nb} &> {log}"


count_variants_docs[
    "Analysis of variant counts"
] = rules.analyze_variant_counts.output.nb


docs["Count variants"] = count_variants_docs
