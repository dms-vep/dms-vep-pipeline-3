"""Functions for ``common.smk``."""


def format_altair_html_chart_params(wc):
    """Get parameters used in ``format_altair_html``."""
    chart = wc.chart
    if m := re.fullmatch(
        "results/func_effects/averages/(?P<condition>.+)_(?P<pheno>func|latent)_effects",
        chart,
    ):
        condition = m.group("condition")
        pheno = m.group("pheno")
        title = (
            func_effects_config["avg_func_effects"][condition]["title"]
            + {
                "func": " (functional score)",
                "latent": " (latent phenotype)",
            }[pheno]
        )
        legend = (func_effects_config["avg_func_effects"][condition]["legend"],)
    elif m := re.fullmatch(
        "results/func_effect_shifts/averages/(?P<comparison>.+)_shifts",
        chart,
    ):
        comparison = m.group("comparison")
        title = avg_func_effect_shifts[comparison]["title"]
        legend = avg_func_effect_shifts[comparison]["legend"]
    elif m := re.fullmatch(
        "results/antibody_escape/averages/(?P<antibody>.+)_mut_(?:escape|icXX)",
        chart,
    ):
        antibody = m.group("antibody")
        title = avg_antibody_escape_config[antibody]["title"]
        legend = avg_antibody_escape_config[antibody]["legend"]
    else:
        raise ValueError(f"could not get formatting params for {chart=}")
    return {"title": title, "legend": legend}
