"""Compare results to expected results (for testing purposes)."""

import os

import pandas as pd


# specify CSVs, merge columns, compare columns
to_compare = [
    (
        "summaries/summary_of_all_phenotypes.csv",
        ["site", "mutant"],
        [
            "monoclonal antibodies escape",
            "spike mediated entry",
            "mock receptor affinity",
        ],
    ),
    (
        "summaries/summary_of_all_phenotypes_per_antibody_escape.csv",
        ["site", "mutant", "antibody"],
        ["escape"],
    ),
    (
        "summaries/summary_of_just_entry_and_binding.csv",
        ["site", "mutant"],
        ["spike mediated entry", "mock receptor affinity"],
    ),
]

for f, merge_cols, compare_cols in to_compare:
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

    # merge, allowing a few rows to be absent in one versus another
    max_diff_rows = 3
    assert len(actual) == len(actual[merge_cols].drop_duplicates())
    assert len(expected) == len(expected[merge_cols].drop_duplicates())
    if abs(len(actual) - len(expected)) > max_diff_rows:
        raise ValueError(f"{len(actual)=}, {len(expected)=}")
    merged = actual[merge_cols + compare_cols].merge(
        expected[merge_cols + compare_cols], on=merge_cols, validate="one_to_one"
    )
    if len(merged) + max_diff_rows < max(len(actual), len(expected)):
        raise ValueError(f"{len(actual)=}, {len(expected)=}, {len(merged)=}")

    for compare_col in compare_cols:
        corr = merged[f"{compare_col}_x"].corr(merged[f"{compare_col}_y"])
        print(f"Correlations for {f}, {compare_col} is {corr:.3f}")
        if corr < 0.99:
            raise ValueError(f"Correlations for {f}, {compare_col} is {corr:.3f}")
