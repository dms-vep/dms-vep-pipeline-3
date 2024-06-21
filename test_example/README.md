# Test example
A small test example. Run with:

    snakemake -j <n_jobs> --software-deployment-method conda -s ../Snakefile

The expected results for some key files (from prior running of test example) are in [./expected_results/](expected_results), and the script [compare_results_to_expected.py](compare_results_to_expected.py) compares the actual results generated in [./results/](results) by running the pipeline to those expected results.
This is designed for testing.
