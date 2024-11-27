"""``snakemake`` files with rules for counting variants from barcode sequencing."""

# Names and values of files to add to docs
count_variants_docs = {
    "Analysis notebooks": collections.defaultdict(dict),
    "Data files": collections.defaultdict(dict),
}


if (
    "use_precomputed_barcode_counts" in config
    and config["use_precomputed_barcode_counts"]
):
    for sample in barcode_runs["sample"]:
        countsfile = f"results/barcode_counts/{sample}_counts.csv"
        if not os.path.isfile(countsfile):
            raise OSError(
                f"`use_precomputed_barcode_counts` is True, but no file {countsfile}"
            )

else:

    rule count_barcodes:
        """Count barcodes for each sample."""
        input:
            fastq_R1=lambda wc: barcode_runs.set_index("sample").at[
                wc.sample, "fastq_R1"
            ],
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
    count_variants_docs["Data files"]["Barcode count CSVs"][sample] = os.path.join(
        "results/barcode_counts",
        f"{sample}_counts.csv",
    )


rule analyze_variant_counts:
    """Analysis of the counts of variants."""
    input:
        expand(
            "results/barcode_counts/{sample}_counts.csv", sample=barcode_runs["sample"]
        ),
        expand(
            "results/barcode_counts/{sample}_invalid.csv",
            sample=barcode_runs["sample"],
        ),
        expand(
            "results/barcode_counts/{sample}_fates.csv", sample=barcode_runs["sample"]
        ),
        gene_sequence_codon=config["gene_sequence_codon"],
        codon_variants=config["codon_variants"],
        site_numbering_map_csv=config["site_numbering_map"],
        barcode_runs_csv=config["barcode_runs"],
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
        """
        papermill {input.nb} {output.nb} \
            -p barcode_runs_csv {input.barcode_runs_csv} \
            -p site_numbering_map_csv {input.site_numbering_map_csv} \
            -p codon_variants {input.codon_variants} \
            -p gene_sequence_codon {input.gene_sequence_codon} \
            &> {log}
        """


if (
    "use_precomputed_barcode_counts" in config
    and not config["use_precomputed_barcode_counts"]
):
    count_variants_docs["Analysis notebooks"][
        "Analysis of variant counts"
    ] = rules.analyze_variant_counts.output.nb


docs["Count barcodes for variants"] = count_variants_docs
