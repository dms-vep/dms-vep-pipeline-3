"""Compare results to expected results (for testing purposes)."""

import os

import pandas as pd


to_compare = [
    "summaries/summary_of_all_phenotypes.csv",
    "summaries/summary_of_all_phenotypes_per_antibody_escape.csv",
    "summaries/summary_of_just_entry_and_binding.csv",
    "summaries/summary_of_just_entry_and_binding_per_antibody_escape.csv",
]

for f in to_compare:
    print(f"\nComparing {f}")
    actual = (
        pd.read_csv(os.path.join("results", f))
        .sort_values(["site", "mutant"])
        .reset_index(drop=True)
    )
    expected = (
        pd.read_csv(os.path.join("expected_results", f))
        .sort_values(["site", "mutant"])
        .reset_index(drop=True)
    )

    print("Comparing actual and expected results...")
    pd.testing.assert_frame_equal(
        actual, expected, check_exact=False, atol=5e-2, rtol=1e-1
    )
    print("Actual and expected results are sufficiently similar.")
