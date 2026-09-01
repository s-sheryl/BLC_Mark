"""Validation of Phase 3 differential-expression inputs for Phase 4."""

from pathlib import Path

import pandas as pd

from .exceptions import EvidenceInputError


REQUIRED_PHASE3_COLUMNS = {
    "gene_id",
    "tested",
    "effect_size",
    "effect_size_label",
    "raw_p_value",
    "adjusted_p_value",
    "significant",
    "missing_reason",
}


def validate_phase3_results(path: str | Path) -> pd.DataFrame:
    """
    Load and validate a Phase 3 differential-expression results file.

    Parameters
    ----------
    path:
        Path to a Phase 3 CSV results file.

    Returns
    -------
    pandas.DataFrame
        Validated Phase 3 results.

    Raises
    ------
    EvidenceInputError
        If the file does not exist, cannot be read, is empty, has missing
        required columns, contains missing gene identifiers, duplicate gene
        identifiers, or untested rows.
    """
    path = Path(path)

    if not path.exists():
        raise EvidenceInputError(f"Phase 3 results file does not exist: {path}")

    if not path.is_file():
        raise EvidenceInputError(f"Phase 3 results path is not a file: {path}")

    try:
        dataframe = pd.read_csv(path)
    except Exception as exc:
        raise EvidenceInputError(
            f"Unable to read Phase 3 results file: {path}"
        ) from exc

    if dataframe.empty:
        raise EvidenceInputError(
            f"Phase 3 results file contains no rows: {path}"
        )

    missing_columns = REQUIRED_PHASE3_COLUMNS - set(dataframe.columns)

    if missing_columns:
        formatted = ", ".join(sorted(missing_columns))
        raise EvidenceInputError(
            f"Phase 3 results file is missing required columns: {formatted}"
        )

    if dataframe["gene_id"].isna().any():
        raise EvidenceInputError(
            "Phase 3 results contain missing gene_id values."
        )

    gene_ids = dataframe["gene_id"].astype(str).str.strip()

    if (gene_ids == "").any():
        raise EvidenceInputError(
            "Phase 3 results contain blank gene_id values."
        )

    if gene_ids.duplicated().any():
        duplicated = sorted(gene_ids[gene_ids.duplicated()].unique())
        raise EvidenceInputError(
            "Phase 3 results contain duplicate gene_id values: "
            + ", ".join(duplicated[:10])
        )

    tested = dataframe["tested"]

    if tested.dtype == bool:
        untested_mask = ~tested
    else:
        untested_mask = (
            tested.astype(str)
            .str.strip()
            .str.lower()
            .ne("true")
        )

    if untested_mask.any():
        raise EvidenceInputError(
            "Phase 3 results contain one or more untested genes."
        )

    dataframe = dataframe.copy()
    dataframe["gene_id"] = gene_ids

    return dataframe