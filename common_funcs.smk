"""Functions for ``common.smk``."""


def format_altair_html_chart_params(wc):
    """Get parameters used in ``format_altair_html``."""
    chart = wc.chart
    if m := re.fullmatch(
        "func_effects/averages/(?P<condition>.+)_(?P<pheno>func|latent)_effects",
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
        "func_effect_shifts/averages/(?P<comparison>.+)_shifts",
        chart,
    ):
        comparison = m.group("comparison")
        title = avg_func_effect_shifts[comparison]["title"]
        legend = avg_func_effect_shifts[comparison]["legend"]
    elif m := re.fullmatch(
        f"(?P<assay>{'|'.join(assays)})/"
        + "averages/(?P<antibody>.+)_mut_(?:effect|icXX)",
        chart,
    ):
        assay = m.group("assay")
        antibody = m.group("antibody")
        title = avg_assay_config[assay][antibody]["title"]
        legend = avg_assay_config[assay][antibody]["legend"]
    elif m != re.fullmatch("summaries/summary_(?:overlaid|faceted)", chart):
        title = summary_config["title"]
        legend = summary_config["legend"]
    else:
        raise ValueError(f"could not get formatting params for {chart=}")
    return {"title": title, "legend": legend}
