"""Compare results to expected results (for testing purposes)."""

import os

import pandas as pd


to_compare = [
    "summaries/summary_of_all_phenotypes.csv",
    "summaries/summary_of_all_phenotypes_per_antibody_escape.csv",
    "summaries/summary_of_just_entry_and_binding.csv",
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

    # check if a substantial fraction of rows are missing, but let there
    # be one or two missing
    if abs(len(actual) - len(expected)) > 3:
        raise ValueError(f"{len(actual)=}, {len(expected)=}")
    else:
        merge_cols = ["site", "mutant"]
        if "antibody" in actual.columns:
            merge_cols.append("antibody")
        assert len(actual) == len(actual[merge_cols].drop_duplicates())
        assert len(expected) == len(expected[merge_cols].drop_duplicates())
        actual = actual.merge(expected[merge_cols], on=merge_cols, validate="1:1")
        expected = expected.merge(actual[merge_cols], on=merge_cols, validate="1:1")

    print("Comparing actual and expected results...")
    pd.testing.assert_frame_equal(
        actual, expected, check_exact=False, atol=0.1, rtol=0.2
    )
    print("Actual and expected results are sufficiently similar.")
