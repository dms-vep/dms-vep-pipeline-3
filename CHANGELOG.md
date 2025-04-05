# CHANGELOG

### 3.22.0
- Allow arbitrary strings as site numbers:
  - Updated `alignparse` to 0.7.0 which enables the use of arbitrary strings as site numbers (eg, `57(E2)` for instance).
  - Made other modifications to code to allow arbitrary numbers. This changes the sorting of the output for some files (eg, the mutation functional effects).

### 3.21.0
- Added `latent_effects_regularization` under `global_epistasis_params` in `func_scores` in `func_effect_configs.yaml`. This is a **backward-incompatible change**, you now **must** specify this parameter. It is used to regularize the mutational effects in the latent space, which leads to better fitting if a very small regularization value is used (see [here](https://github.com/matsengrp/multidms/issues/168) and [here](https://github.com/dms-vep/dms-vep-pipeline-3/issues/182)). If you set this parameter to zero, the results will be identical to earlier versions of this pipeline. If you set to a very small value like 1e-7 (as recommended) it will lead to slightly different results than earlier versions of this pipeline. In some cases, you could want regularization even more than 1e-7 (see [here](https://github.com/dms-vep/dms-vep-pipeline-3/issues/182#issuecomment-2715852069).

#### version 3.20.1
- Fix [minor issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/176) w embedding Altair plots in VitePress homepage.

### version 3.20.0
- Various updates to the `summary` plots:
  + Add option `no_mean_lineplot` to not show antibody-escape mean lineplot (see [this issue](https://github.com/dms-vep/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/issues/136)).
  + Add option `lineplot_antibody_label_loc` ("right" or "top") for labels on antibody-escape lineplots ((see [this issue](https://github.com/dms-vep/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS/issues/136)).
  + Add site labels to lineplots when only one antibody per group (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/177))
  + Add `scale_lineplot_height` option to adjust height of lineplots.
  + Add `selectable_per_antibody_heatmap` option to make antibody escape heatmap selectable for which antibody is shown.

#### version 3.19.3
- Fix bug in `func_effect_diffs` tooltips introduced in 3.19.2.

#### version 3.19.2
- Make `func_effect_diffs` work when there are many selections by making correlation heatmap rather than many scatters, and showing tooltips as list of per-selection values. Also makes `per_selection_tooltips` not needed in config for `func_effect_diffs`.

#### version 3.19.1
- Fixed [bug](https://github.com/dms-vep/dms-vep-pipeline-3/issues/163) introduced in version 3.19.0 that causes pip installation of `alignparse==0.6.3` to fail.

### version 3.19.0
- Added the ability to average results (functional effects, functional effect shifts, antibody escape) for different experiments keeping only different regions (sites) for each selection. This is useful if the experiments are done on different halves of a gene. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/160).
- The `func_effects_global_epistasis` notebooks now make plots showing correlation between mutation effects inferred for mutations by the global epistasis model versus the average effects of these mutations among all single-mutant variants. This helps with evaluating questions like [this](https://github.com/dms-vep/dms-vep-pipeline-3/issues/158) about whether single mutant estimates are being affected by the global epistasis model.
- Increased the number of iterations and the stringency of the convergence criteria for the `multidms` global epistasis model fitting per the suggestions [here](https://github.com/dms-vep/dms-vep-pipeline-3/issues/158#issuecomment-2499071819).
- Update `conda` environment to latest versions of `altair` (5.4 -> 5.5), `baltic` (0.2.2 -> 0.3), `biopython` (1.83 -> 1.84), `entrez-direct` (16.2 -> 22.4), `markdown` (3.4 -> 3.6), `matplotlib` (3.8 -> 3.9). `minimap2` (2.26 -> 2.28), `plotnine` (0.12 -> 0.14), `python` (3.11 -> 3.12)`, `ruamel.yaml` (0.17 -> 0.18), `snakemake` (8.20.6 -> 8.25.3), `alignparse` (0.6.2 -> 0.6.3), `polyclonal` (6.12 -> 6.14), and `multidms` (0.3.3 -> 0.4.2). Note that `multidms` 0.4.2 fixes a [this issue](https://github.com/matsengrp/multidms/issues/165) where `altair` was being downgraded to 5.1.2. Note that `polyclonal` 6.14 fixes [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/107) where the `minimum max <stat> of site` slider did not work properly with hidden sites.
- Made some internal changes to code to adjust to version update of `ruamel.yaml`, as new version drops `round_trip_dump` and `safe_load`.
- Some notebooks previously read input parameters from the YAML configuration file, which is non-ideal as it meant that `snakemake` was not explicitly tracking what was passed to them. Fixed this: now all parameters they use are passed via `papermill` parameterization.
- Fixed [bug](https://github.com/dms-vep/dms-vep-pipeline-3/issues/164) introduced in version 3.14 that causes `analyze_variant_counts` notebook to be excluded from docs in some cases.
- Update node version for homepage deployment and add note to README about how deployment version must match that used to create `package-lock.json` (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/165)).

### version 3.18.0
- Add `no_heatmap` option to `summaries_config` to not show some heatmaps. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/155).

### version 3.17.0
- Allow JSON (`*.json`) files to be listed in the auto-rendered DOCS in the same way CSV and FASTA files can. Useful for linking to `dms-viz` JSONs.

#### version 3.16.3
- Update `snakemake`, `altair`, and `upsetplot` versions in `environment.yml`.

#### version 3.16.2
- Remove `defaults` channel from `environment.yml`. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/151).

#### version 3.16.1
- Update to `polyclonal` 6.12 in order to work with latest `binarymap`; otherwise the update to the newest `binarymap` version caused error.

### version 3.16.0
- Add `duplicate_fastq_R1` flag to `config.yaml` to allow the user to specify what happens if a FASTQ is duplicated among samples. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/147).

### version 3.15.0
- Update some packages in `conda` env (fixes [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/144)):
  - `pandas` to 2.2
  - `multidms` to 0.3.3
  - `altair` to 5.3
  - `snakemake` to 8.14
  - `neutcurve` to 2.1
  - `papermill` to 2.6
- Added test that compares actual results to expected ones for `test_example` better testing of pipeline updates.

#### version 3.14.1
- Make change in version 3.14.0 backward compatible by assuming `use_precomputed_barcode_counts` is `False` if not specified.

### version 3.14.0
- Add option to use pre-computed barcode counts (`use_precomputed_barcode_counts`) so pipeline can be entirely run from pre-calculated barcode-variant table and barcode counts within repo. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/141).

### version 3.13.0
- Add replicate scatter plots in the notebooks that average mutation effects. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/137).

### version 3.12.0
- Make nicer version of the scatter plot produced by `func_effect_diffs` that compares the functional effects in the two conditions, and save it to the default docs.

### version 3.11.0
- Allow more customization of summary plots, such as by specifying CSVs explicitly. This is a **backward-incompatible** change in how specify the YAML configuration for the summary plots, in that now the CSV file is specified and all non-escape phenotypes are grouped together rather than as separate keys.
- Fix bug in summary plot line plot zooming. Now sites in line plot are always aligned with overlay bar prior to zooming, and missing sites are shown empty. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/124).
- Fix bug in setting of limits in summary plots with `min_at_least`.

### version 3.10.0
- Remove the titles and legends from interactive figures (see [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/121)). These were not really being used, and with the new VitePress option that is a better way to provide detailed information around figures. This change:
  - removes the `title` and `legend` keys associated with various figures and plots in the configuration YAMLs.
  - the pipeline no longer makes the `<path>_nolegend.html` versions of files as the final version (`<path>.html`) now does not have a legend, so corresponding remove the `format_altair_html` rule. This removes the following files:
    - `common.smk`
    - `common_funcs.smk`
    - `scripts/format_altair_html.py`
  - Remove `bs4` from `environment.yml` as it is no longer needed.

### version 3.9.0
- Allow `mutation_annotations` file that provides annotations for specific mutations, such as how many nucleotide mutations are required to generate them (see [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/105)). In the test example, this is implemented by specifying `nt changes to codon`. These annotations can then be used as filters by specifying the indicated columns in the configuration for the average across replicates plots for antibody escape and functional effects, as well as in the summaries.
- Bug fix to alignment of line plot and scale bar and some filtering bugs in summary plots

#### version 3.8.1
- Fix change to add VitePress homepage backward-compatible as originally intended by not building homepage (versus raising error) in VitePress related data totally ommitted from config.

### version 3.8.0
- Add code and instructions on how to optionally build a nicer and more customized documentation using VitePress for hosting on GitHub Pages.
- Update the version of `mafft` from `7.520` to `7.525` to avoid channel priority issues.

### version 3.7.0
- Allow multiple summaries to be specified in `summaries_config.yaml`. Addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/110). Note that the name of the output files in `./results/summaries` has slightly changed, and as a result you will need to update the `.gitignore` slightly. When migrating to this new version, full delete your old `./results/summaries` subdirectory to remove obsolete naming, then re-run. Also, the table of contents in the output HTML for GitHub pages now labels the summaries slightly differently.
- Lint Jupyter notebooks by fixing `ruff.toml` (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/111)).

#### version 3.6.3
- Enable `func_effect_diffs` heat map to be filtered by mutation effects on each individual type of functional effect (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/108)). This change actually just updates the `func_effects_config.yml` for the test example, not the actual code.

#### version 3.6.2
- Pipeline can run correctly when `antibody_escape_config` is null (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/102)).

#### version 3.6.1
- Better error messages with problems with `barcode_runs`

### version 3.6.0
- Add `func_effect_diffs` option to compare differences among functional effects.
- Update environment:
  - update `polyclonal` to 6.11 (addresses [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/96))
  - update `neutcurve` to 1.1.2
  - update `altair` to 5.2.0
  - update `biopython` to 1.83
  - update `pandas` to 2.1 and add `pyarrow`. Did not update to `pandas` 2.2 due to [this issue](https://github.com/matsengrp/multidms/issues/128).
  - update to `seaborn` 0.13
  - update to `snakemake` 8.3. **Note that this means the recommended usage now changes from `--use-conda` to `--software-deployment-method conda`.**
- sort rows in prob escape values for consistent output, may very slightly change some of the fit antibody-escape values

#### version 3.5.6
- Fit each lasso weight independently in `multidms` comparisons (see [this issue](https://github.com/dms-vep/dms-vep-pipeline-3/issues/92)).

#### version 3.5.5
- Fix `PeriodicWildCardError` that was raised by `snakemake` for the `_nolegend` output HTML plots for certain wildcards.
- Fix how notebook linting done with `ruff` in tests.

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
