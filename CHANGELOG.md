# CHANGELOG

#### version 3.2.3
- Add `other_target_files` variable: files can be added to this (eg, in `custom_rules.smk`) and then the pipeline will ensure these files are also created. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/39).

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
