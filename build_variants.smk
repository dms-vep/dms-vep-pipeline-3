"""``snakemake`` files with rules for building variants."""

build_variants_docs = {}


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
            "results/logs/get_prebuilt_variants.txt",
        shell:
            """
            curl -o {output.variants} {params.variants_url} &> {log}
            curl -o {output.geneseq} {params.geneseq_url} &>> {log}
            """

else:
    pacbio_runs = pd.read_csv(config["pacbio_runs"])
    pacbio_runs_required_columns = {"library", "run", "fastq"}
    if not pacbio_runs_required_columns.issubset(pacbio_runs.columns):
        raise ValueError(
            f"{pacbio_runs.columns=} lacks {pacbio_runs_required_columns=}"
        )
    if len(pacbio_runs) != pacbio_runs["run"].nunique():
        raise ValueError(
            f"duplicates in 'run' column of `pacbio_runs`\n\n{pacbio_runs}"
        )

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
            "results/logs/align_parse_PacBio_ccs_{pacbioRun}.txt",
        script:
            "scripts/align_parse_PacBio_ccs.py"

    rule analyze_pacbio_ccs:
        """Analyze PacBio CCSs and get ones that align to amplicons of interest."""
        input:
            expand(
                rules.align_parse_PacBio_ccs.output.outdir,
                pacbioRun=pacbio_runs["run"],
            ),
            pacbio_amplicon=config["pacbio_amplicon"],
            pacbio_amplicon_specs=config["pacbio_amplicon_specs"],
            pacbio_runs_csv=config["pacbio_runs"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/analyze_pacbio_ccs.ipynb"
            ),
        output:
            csv="results/process_ccs/CCSs_aligned_to_amplicon.csv",
            nb="results/notebooks/analyze_pacbio_ccs.ipynb",
        conda:
            "environment.yml"
        log:
            "results/logs/analyze_pacbio_ccs.txt",
        shell:
            """
            papermill {input.nb} {output.nb} \
                -p pacbio_amplicon {input.pacbio_amplicon} \
                -p pacbio_amplicon_specs {input.pacbio_amplicon_specs} \
                -p pacbio_runs_csv {input.pacbio_runs_csv} \
                &> {log}
            """

    rule build_pacbio_consensus:
        """Build PacBio consensus sequences for barcodes."""
        input:
            rules.analyze_pacbio_ccs.output.csv,
            gene_sequence_codon=config["gene_sequence_codon"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/build_pacbio_consensus.ipynb"
            ),
        output:
            nt_variants="results/variants/nt_variants.csv",
            nb="results/notebooks/build_pacbio_consensus.ipynb",
        params:
            max_error_rate=config["max_ccs_error_rate"],
            params_yaml=yaml_str(
                {
                    "consensus_params": config["consensus_params"],
                    "variant_tags": config["variant_tags"],
                }
            ),
        conda:
            "environment.yml"
        log:
            "results/logs/build_pacbio_consensus.txt",
        shell:
            """
            papermill {input.nb} {output.nb} \
                -p gene_sequence_codon {input.gene_sequence_codon} \
                -p max_error_rate {params.max_error_rate} \
                -y "{params.params_yaml}" \
                &> {log}
            """

    rule build_codon_variants:
        """Build codon-variant table."""
        input:
            rules.build_pacbio_consensus.output.nt_variants,
            gene_sequence_codon=config["gene_sequence_codon"],
            gene_sequence_protein=config["gene_sequence_protein"],
            site_numbering_map_csv=config["site_numbering_map"],
            neut_standard_barcodes=config["neut_standard_barcodes"],
            mutation_design_classification_csv=config["mutation_design_classification"]["csv"],
            nb=os.path.join(
                config["pipeline_path"], "notebooks/build_codon_variants.ipynb"
            ),
        output:
            codon_variants=config["codon_variants"],
            nb="results/notebooks/build_codon_variants.ipynb",
        params:
            mutation_design_classification_site_col=config["mutation_design_classification"]["site_col"],
        conda:
            "environment.yml"
        log:
            "results/logs/build_codon_variants.txt",
        shell:
            """
            papermill {input.nb} {output.nb} \
                -p gene_sequence_codon {input.gene_sequence_codon} \
                -p gene_sequence_protein {input.gene_sequence_protein} \
                -p site_numbering_map_csv {input.site_numbering_map_csv} \
                -p neut_standard_barcodes {input.neut_standard_barcodes} \
                -p mutation_design_classification_csv {input.mutation_design_classification_csv} \
                -p mutation_design_classification_site_col {params.mutation_design_classification_site_col} \
                -p codon_variants {output.codon_variants} \
                &> {log}
            """

    build_variants_docs["Analysis notebooks"] = {
        "Analysis of PacBio CCSs": rules.analyze_pacbio_ccs.output.nb,
        "Building barcoded variant consensus sequences": rules.build_pacbio_consensus.output.nb,
        "Building of codon-variant table": rules.build_codon_variants.output.nb,
    }


# Names and values of files to add to docs
build_variants_docs["Data files"] = {
    "parental codon sequence": config["gene_sequence_codon"],
    "parental protein sequence": config["gene_sequence_protein"],
    "barcode to codon-variant lookup table": config["codon_variants"],
}

docs["Barcode to codon-variant lookup table"] = build_variants_docs
