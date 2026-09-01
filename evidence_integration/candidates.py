"""Candidate extraction for BLC Mark Phase 4 evidence integration."""

from pathlib import Path

import pandas as pd

from .validation import validate_phase3_results


def extract_significant_candidates(path: str | Path) -> pd.DataFrame:
    """
    Extract significant Phase 3 genes for Phase 4 evidence integration.

    Parameters
    ----------
    path:
        Path to a validated Phase 3 differential-expression CSV file.

    Returns
    -------
    pandas.DataFrame
        Significant candidate genes only.

    Notes
    -----
    This function does not rank candidates and does not introduce any
    additional effect-size threshold. It preserves the Phase 3 significance
    decision exactly as recorded in the input results.
    """
    dataframe = validate_phase3_results(path)

    significant = dataframe["significant"]

    if significant.dtype == bool:
        mask = significant
    else:
        mask = (
            significant.astype(str)
            .str.strip()
            .str.lower()
            .eq("true")
        )

    candidates = dataframe.loc[mask].copy()

    candidates = candidates.sort_values(
        by=["adjusted_p_value", "gene_id"],
        ascending=[True, True],
        na_position="last",
        kind="stable",
    ).reset_index(drop=True)

    return candidates