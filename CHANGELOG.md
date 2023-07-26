# CHANGELOG

### version 3.1.0
- Update `multidms` to version 0.2.0

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
