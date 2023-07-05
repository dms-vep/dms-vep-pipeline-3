"""``snakemake`` files with rules for calculating functional effects."""


# read the config for func effects
with open(config["func_effects_config"]) as f:
    func_effects_config = yaml.safe_load(f)

# get configuration for functional scores and make sure all samples defined
func_scores = func_effects_config["func_scores"]
for selection_name, selection_d in func_scores.items():
    for s in ["post_selection_sample", "pre_selection_sample"]:
        if selection_d[s] not in sample_to_library:
            raise ValueError(
                f"{s} for {selection_name} of {selection_d[s]} not in barcode_runs"
            )

# Names and values of files to add to docs
func_effects_docs = collections.defaultdict(dict)


rule func_scores:
    """Compute functional scores for variants."""
    input:
        unpack(
            lambda wc: {
                s: f"results/barcode_counts/{func_scores[wc.selection][s]}_counts.csv"
                for s in ["post_selection_sample", "pre_selection_sample"]
            }
        ),
        codon_variants=config["codon_variants"],
        gene_sequence_codon=config["gene_sequence_codon"],
        site_numbering_map=config["site_numbering_map"],
    output:
        csv="results/func_scores/{selection}_func_scores.csv.gz",
    params:
        func_score_params=lambda wc: func_scores[wc.selection]["func_score_params"],
        # script will throw error if pre_library and post_library differ
        libraries=lambda wc: {
            s: sample_to_library[func_scores[wc.selection][s]]
            for s in ["post_selection_sample", "pre_selection_sample"]
        },
    conda:
        "environment.yml"
    log:
        "results/logs/func_scores_{selection}.txt",
    script:
        "scripts/func_scores.py"


for s in func_scores:
    func_effects_docs["Functional scores"][s] = f"results/func_scores/{s}_func_scores.csv.gz"


docs["Functional effects of mutations"] = func_effects_docs
