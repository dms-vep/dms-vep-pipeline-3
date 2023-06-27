"""``snakemake`` files with rules for calculating functional effects."""


# read the config for func effects
with open(config["func_effects_config"]) as f:
    func_effects_config = yaml.safe_load(f)

# Names and values of files to add to docs
func_effects_docs = {}


rule func_scores:
    """Compute functional scores for variants."""
    input:
        unpack(
            lambda wc: {
                s: (
                    "results/barcode_counts/"
                    + func_effects_config["func_effects_selections"][wc.condition][
                        "selections"
                    ][wc.selection][s]
                    + "_counts.csv"
                )
                for s in ["post_selection_sample", "pre_selection_sample"]
            }
        ),
    output:
        csv="results/func_scores/{condition}/{selection}.csv",
    params:
        lambda wc: func_effects_config["func_effects_selections"][wc.condition][
            "selections"
        ][wc.selection],
    conda:
        "environment.yml"
    log:
        "results/logs/func_scores_{condition}_{selection}.txt",
    script:
        "scripts/func_scores.py"


for condition, condition_d in func_effects_config["func_effects_selections"].items():
    for selection in condition_d["selections"]:
        func_effects_docs[
            f"{condition} {selection}"
        ] = f"results/func_scores/{condition}/{selection}.csv"


docs["Functional effects of mutations"] = func_effects_docs
