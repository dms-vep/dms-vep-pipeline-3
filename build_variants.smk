"""``snakemake`` files with rules for building variants."""


# Names and values of files to add to docs
build_variants_docs = {
    "parental codon sequence": config["gene_sequence_codon"],
    "parental protein sequence": config["gene_sequence_protein"],
    "barcode to codon-variant lookup table": config["codon_variants"],
}
    
rule translate_geneseq:
    """Translate gene sequence into protein sequence."""
    input:
        gene=config["gene_sequence_codon"],
    output:
        prot=config["gene_sequence_protein"],
    conda:
        "environment.yml"
    log:
        "results/logs/translate_geneseq.txt",
    script:
        "scripts/translate_geneseq.py"


if (config["prebuilt_variants"] is None) != (config["prebuilt_geneseq"] is None):
    raise ValueError("specify both or neither `prebuilt_geneseq` / `prebuilt_variants`")

elif config["prebuilt_variants"] and config["prebuilt_geneseq"]:

    rule get_prebuilt_variants:
        """Get pre-built variants and gene sequence."""
        output:
            variants=config["codon_variants"],
            geneseq=config["gene_sequence_codon"],
        params:
            variants_url=config["prebuilt_variants"],
            geneseq_url=config["prebuilt_geneseq"],
        conda:
            "environment.yml"
        log:
            "results/logs/get_prebuilt_variants.txt"
        shell:
            """
            curl -o {output.variants} {params.variants_url} &> {log}
            curl -o {output.geneseq} {params.geneseq_url} &>> {log}
            """

else:

    pacbio_runs = pd.read_csv(config["pacbio_runs"])
    pacbio_runs_required_columns = {"library", "run", "fastq"}
    if not pacbio_runs_required_columns.issubset(pacbio_runs.columns):
        raise ValueError(f"{pacbio_runs.columns=} lacks {pacbio_runs_required_columns=}")
    if len(pacbio_runs) != pacbio_runs["run"].nunique():
        raise ValueError(f"duplicates in 'run' column of `pacbio_runs`\n\n{pacbio_runs}")

    rule gene_sequence:
        """Get sequence of gene from PacBio amplicon."""
        input:
            gb=config["pacbio_amplicon"],
        output:
            codon=config["gene_sequence_codon"],
        conda:
            "environment.yml"
        log:
           "results/logs/gene_sequence.txt",
        script:
            "scripts/gene_sequence.py"


    rule align_parse_PacBio_ccs:
        """Align and parse PacBio CCS FASTQ file."""
        input:
            fastq=lambda wc: pacbio_runs.set_index("run").at[wc.pacbioRun, "fastq"],
            amplicon=config["pacbio_amplicon"],
            specs=config["pacbio_amplicon_specs"],
        output:
            outdir=directory("results/process_ccs/{pacbioRun}"),
        conda:
            "environment.yml"
        log:
            "results/logs/align_parse_PacBio_ccs_{pacbioRun}.txt"
        script:
            "scripts/align_parse_PacBio_ccs.py"


    rule analyze_pacbio_ccs:
        """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
        input:
            expand(
                rules.align_parse_PacBio_ccs.output.outdir,
                pacbioRun=pacbio_runs["run"],
            ),
            config["pacbio_amplicon"],
            config["pacbio_amplicon_specs"],
            nb=os.path.join(config["pipeline_path"], "notebooks/analyze_pacbio_ccs.ipynb"),
        output:
            csv="results/process_ccs/CCSs_aligned_to_amplicon.csv",
            nb="results/notebooks/analyze_pacbio_ccs.ipynb",
        conda:
            "environment.yml"
        log:
            "results/logs/analyze_pacbio_ccs.txt"
        shell:
            "papermill {input.nb} {output.nb} &> {log}"

    build_variants_docs["Analysis of PacBio CCSs"] = rules.analyze_pacbio_ccs.output.nb


    rule build_pacbio_consensus:
        """Build PacBio consensus sequences for barcodes."""
        input:
            rules.analyze_pacbio_ccs.output.csv,
            config["gene_sequence_codon"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/build_pacbio_consensus.ipynb"
            ),
        output:
            nt_variants="results/variants/nt_variants.csv",
            nb="results/notebooks/build_pacbio_consensus.ipynb",
        params:
            config["max_ccs_error_rate"],
            config["consensus_params"],
        conda:
            "environment.yml"
        log:
            "results/logs/build_pacbio_consensus.txt"
        shell:
            "papermill {input.nb} {output.nb} &> {log}"

    build_variants_docs[
        "Building of consensus sequences for barcoded variants"
    ] = rules.build_pacbio_consensus.output.nb


    rule build_codon_variants:
        """Build codon-variant table."""
        input:
            rules.build_pacbio_consensus.output.nt_variants,
            config["gene_sequence_codon"],
            config["gene_sequence_protein"],
            config["site_numbering_map"],
            config["mutation_design_classification"]["csv"],
            config["neut_standard_barcodes"],
            nb=os.path.join(config["pipeline_path"], "notebooks/build_codon_variants.ipynb"),
        output:
            config["codon_variants"],
            nb="results/notebooks/build_codon_variants.ipynb",
        params:
            config["mutation_design_classification"],
        conda:
            "environment.yml"
        log:
            "results/logs/build_codon_variants.txt"
        shell:
            "papermill {input.nb} {output.nb} &> {log}"

    build_variants_docs[
        "Building of codon-variant table"
    ] = rules.build_codon_variants.output.nb

docs["Barcode to codon-variant lookup table"] = build_variants_docs
