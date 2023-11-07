# CHANGELOG

#### version 3.5.4
- Fix bug in reporting pre-selection counts cutoff in `func_scores`. See [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/87). The cutoff was applied correctly before and that has not changed, this just fixes reporting of it in plots in `analyze_func_scores`.

#### version 3.5.3
- Fix bug in `per_antibody_escape.csv` production in `summary` introduced in version 3.5.1.

#### version 3.5.2
- Fix bug introduced in `summary` when no antibodies in version 3.5.1

#### version 3.5.1
- Fix tooltips in std chart in averaging notebooks when column has a `.` in the name (see [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/81)).
- `summary` rule now creates a per-antibody escape CSV file.

### version 3.5.0
- Updated to `polyclonal` 6.9
- Add options to filter average measurements based on excessive variability in measurements across replicates (see [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/78)):
  - `func_effects_config.yaml` now allows you to specify `effect_std: <value>` under `plot_kwargs: addtl_slider_stats`, remember to also add `addtl_slider_stats_as_max: [effect_std]`. You can also specify `floor_for_effect_std` under `avg_func_effects` to floor before computing the standard deviation.
  - `avg_func_effects.ipynb` now makes plots with `effect_std` slider and also plots distribution of effect standard deviations to help you choose a good initial value for this filter.
  - the results files with the average functional effects now include the `effect_std` column
  - `antibody_escape_config.yaml` now allows you to specify `escape_std: <value>` under `plot_kwargs: addtl_slider_stats`, remember to also add `addtl_slider_stats_as_max: [escape_std]`.
  - `avg_antibody_escape.ipynb` now makes plots with `escape_std` slider and also plots distribution of escape standard deviations to help you choose a good initial value for this filter.
  - the results files with the average antibody escape now include the `escape_std` column
  - `summaries_config.yml` now provides a `le_filter` option which can be used to set a filter on the standard deviations.
- Customize scales in probability escape plots by adding `concentration_scale` and `concentration_title` to assay configuration in `antibody_escape_config.yaml`. See [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/72).

#### version 3.4.10
- Added `baltic` to the environment.

#### version 3.4.9
- Update to `dmslogo` 0.7.0, which fixes [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/74)

#### version 3.4.8
- Update versions of some software in `conda` env:
  - `altair` -> 5.1.2
  - `matplotlib` -> 3.8
  - `papermill` -> 2.4
  - `snakemake` -> 7.32
  - `dmslogo` -> added version 0.6.3
  - `neutcurve` -> added version 0.5.7

#### version 3.4.7
- In `summary` plots you can now specify **either** `max_at_least` or `fixed_max` and `min_at_least` or `fixed_min` for heatmaps.

#### version 3.4.6
- Order antibodies in `summary` plots as ordered in config. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/66).

#### version 3.4.5
- Fix bug that occurs when `barcode_runs` is set to `null` in the `config.yaml`. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/64).

#### version 3.4.4
- Fix bug in `summary` when no `other_assays`; also fix filtering on heatmaps sliders when more than two other properties. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/62).

#### version 3.4.3
- Antibody names passed to `summary` cannot have non-alphanumeric characters. Before this led to a silent and hard to diagnose failure. Now there is an explicit check that throws an interpretable error if any antibodies are assigned for `summary` non-alphanumeric names (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/58)).

#### version 3.4.2
- Fix bug so `summary` works even if no antibody escape data.
- Include legend in `summary` plot when no antibody escape data.

#### version 3.4.1
- Fix bug in site sorting in line plots in `summary.ipynb`

### version 3.4.0
- Add new rules to create summaries of all of the data across assays. This adds the new rules in `summaries.smk` and the new configuration in `summaries_config.yml`.
- Add `show_icXX_in_docs` key to `antibody_escape_config.yml` under `avg_assay`, and then do **not** show ICXX values in final docs if this is `False` or missing. The reason is that these values seem to be highly correlated with the escape values themselves, so there is often not much advantage showing them additionally. Furthermore, they sometime seem to be linearly correlated with validation assays but not with slope 1, making the unit-less escape value perhaps better (and so giving a reason not to show ICXX). This change is **backward incompatible** with respect to how the final HTML docs look: unless you add `show_icXX_in_docs: true` for a given assay average, it will not be shown in docs. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/52).
- Fix bug introduced in version 3.3.0 ([here](https://github.com/dms-vep/dms-vep-pipeline-3/pull/50)) where the information written to the `avg_escape` CSVs was not correct for non-antibody and non-receptor-affinity selections.

#### version 3.3.1
- Improve visual appearance of some site plots by making x-axis labels less crowded.
- In `avg_escape`, plot the correlations using the same filters and metric that is being used for the final average plot. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/46).

### version 3.3.0
- Restructure to allow arbitrary assays in `antibody_escape_config.yml`, not just antibody escape and receptor affinity. **Backwards incompatible**: This restructuring requires you to add a new `assays` key to `antibody_escape_config.yml` that defines the assays being used.

#### version 3.2.5
- Do not show assays in the docs if no samples for that assay.

#### version 3.2.4
- Allow additional site numbering schemes to be defined in addition to `reference_site` and `sequential_site` by retaining as a tooltip any column in `site_numbering_map` that ends in `site`. This change also means you no longer need to specify `sequential_site` under `addtl_tooltip_stats` in `antibody_escape_config.yml`.

#### version 3.2.3
- De-clutter the HTML docs by creating subsections that separate the "Final summary plots", "Analysis notebooks", and "Data files". Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/41).
- Add `other_target_files` variable: files can be added to this (eg, in `custom_rules.smk`) and then the pipeline will ensure these files are also created. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/39).
- Updated to `polyclonal` 6.7, which addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/40) and [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/43).

#### version 3.2.2
- Add `min_filters` to `plot_hide_stats` in `antibody_escape_config.yaml` to only hide/filter robustly measured values.

#### version 3.2.1
- Constrain wildcards in `format_altair_html` to better enable custom rules

### version 3.2.0
- Generalize `antibody_escape` to also allow receptor affinity selections with soluble receptor.
- Changed output of `avg_escape` rules for antibodies to be `mut_effect` rather than `mut_escape`. This is because `effect` can generally refer to other phenotypes too.
- For the average mutation effects and ICXX values, also report the per-model values.
- Upgrade to `polyclonal` 6.6, which now shows (hidden) mutations w deleterious functional effects of other hiding metrics.
- Collapse lists in docs index even if just one entry
- Apply `times_seen` filter when computing correlations in `avg_escape`.
- Internally restructured to re-use a single rule for formatting altair plots by creating `common.smk` and then putting the `format_altair_html` rule there and removing plot-specific versions of that rule.

#### version 3.1.3
- Fixed a bug in `avg_func_effect_shifts` when lasso penalty is in scientific notation.

#### version 3.1.2
- Upgrade to `polyclonal` 6.5.

#### version 3.1.1
- For average plots, change default to only show mutations if present in one more than half of all models / selections / comparisons. This reduces showing of mutations found in a minority of replicates.

### version 3.1.0
- Compute and analyze shifts in functional effects between different conditions computed using `multidms`. This adds the `func_effect_shifts`, `avg_func_effect_shifts`, and `format_avg_func_effect_shifts_chart` rules.
- Update `multidms` to version 0.2.1
- Add `ipywidgets` and `seaborn` to environment
- Improvements to GitHub actions testing: save log files on failure, and don't build `conda` environment twice.

#### version 3.0.8
- Allow underscores in antibody selection names. (Note this changes paths for prob escape files).

#### version 3.0.7
- Add tooltips to plot in `avg_antibody_escape`
- Use log-scale in neutralization vs concentration plots for concentration.

#### version 3.0.6
- Added plot of diversity indices in `analyze_variant_counts`

#### version 3.0.5
- Fix x-axis range of line plots in `avg_antibody_escape`
- Plot x-axis using reference rather than sequential site in `build_codon_variants` and `analyze_variant_counts`

#### version 3.0.4
- Allow the parent gene sequence to not start with ATG (start codon) if it fully translates otherwise.

#### version 3.0.3
- Update to `polyclonal` 6.4

#### version 3.0.2
- Fix error causing crashing of `build_docs` due to leftover line from earlier debugging.
- Update to `polyclonal` 6.3
- Updated to `multidms` 0.1.9

#### version 3.0.1
- Upgrade `multidms` to version 0.1.6, which fixes error with letter-suffixed site numbers.

## version 3.0.0
This is a total re-write of the pipeline in a new repository, replacing versions 1 and 2 that are still available at [https://github.com/dms-vep/dms-vep-pipeline](https://github.com/dms-vep/dms-vep-pipeline).
There are numerous changes to the pipeline structure in version 3 relative to versions 1 and 2, although much of the logic of the pipeline is the same.
