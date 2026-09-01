"""
Purpose:
    Perform structural validation of DE inputs (expression matrix and
    sample metadata) before any statistical computation happens, and
    reconcile expression-matrix samples against metadata samples.

Responsibilities:
    - Load and structurally validate the expression matrix: file
      existence/readability, gene-identifier column presence, sample
      columns present, unique gene identifiers, unique sample column
      names, numeric expression values.
    - Load and structurally validate sample metadata: file
      existence/readability, sample-identifier column presence, group
      column presence, unique sample identifiers, non-empty group
      labels.
    - Match expression-matrix sample columns against metadata sample
      identifiers, reporting matched, expression-only, and
      metadata-only samples without silently discarding any of them.
    - Fail explicitly (raising the appropriate
      differential_expression.exceptions subclass) whenever a
      required condition is not met.

Scope:
    This module performs structural validation only. It does not
    perform statistical testing, does not decide which group is the
    reference or comparison group (that belongs to comparison.py),
    and does not read raw, unprocessed downloaded files -- the
    expression matrix it validates must already be a processed
    matrix (specification Section 3.1).

    This module intentionally implements its own lightweight
    delimiter-aware table loading rather than reusing
    preprocessing_manager.PreprocessingManager.load_dataset(), because
    that method is an instance method of a class that requires a full
    ProjectConfig and is designed around the raw-file structural
    cleanup responsibilities of the preprocessing stage, not around
    validating an already-processed matrix for DE-specific semantics
    (gene-ID/sample-column/group-label checks). The delimiter
    convention (".csv" -> comma, ".tsv"/".txt" -> tab, ".gz" ->
    compressed inner delimiter) matches preprocessing_manager.py's own
    convention for consistency, without importing or duplicating its
    class internals.
"""

from pathlib import Path
from typing import Type

import pandas as pd

from src.differential_expression.exceptions import (
    InvalidExpressionMatrixError,
    InvalidMetadataError,
    SampleMismatchError,
)
from src.differential_expression.models import (
    SampleMatchResult,
    ValidatedExpressionMatrix,
    ValidatedMetadata,
)

VALIDATION_VERSION = "1.0"

_DELIMITER_BY_EXTENSION = {
    ".csv": ",",
    ".tsv": "\t",
    ".txt": "\t",
}

__all__ = [
    "VALIDATION_VERSION",
    "validate_expression_matrix",
    "validate_metadata",
    "match_samples",
]


def _resolve_delimiter(path: Path, error_cls: Type[Exception]) -> str:
    """Infer a delimiter from a file's extension.

    Args:
        path: The file whose delimiter should be inferred.
        error_cls: Exception type to raise if the extension is not
            recognized.

    Returns:
        The delimiter character to use with pandas.read_csv.

    Raises:
        error_cls: If the extension (or, for ".gz", the inner
            extension) is not one of ".csv", ".tsv", ".txt".
    """
    suffixes = path.suffixes
    if path.suffix == ".gz" and len(suffixes) > 1:
        inner_suffix = suffixes[-2]
    else:
        inner_suffix = path.suffix

    if inner_suffix not in _DELIMITER_BY_EXTENSION:
        raise error_cls(
            f"Cannot infer a delimiter for {path}; recognized "
            f"extensions are {sorted(_DELIMITER_BY_EXTENSION)} "
            "(optionally gzip-compressed)."
        )

    return _DELIMITER_BY_EXTENSION[inner_suffix]


def _load_table(
    path: Path,
    context: str,
    error_cls: Type[Exception],
) -> pd.DataFrame:
    """Load a delimited table, raising a DE-specific exception on
    every failure mode instead of letting a generic pandas/OS
    exception propagate.

    Args:
        path: Path to the file to load.
        context: Human-readable label ("Expression matrix" or
            "Metadata") used in error messages.
        error_cls: The DE exception type to raise on failure.

    Returns:
        The loaded table, exactly as read (no cleaning applied).

    Raises:
        error_cls: If the path does not exist, is not a file, cannot
            be parsed, or parses to an empty table.
    """
    if not isinstance(path, Path):
        raise TypeError(f"'{context}' path must be a pathlib.Path, got {type(path).__name__}.")

    if not path.exists():
        raise error_cls(f"{context} file does not exist: {path}")

    if not path.is_file():
        raise error_cls(f"{context} path exists but is not a file: {path}")

    separator = _resolve_delimiter(path, error_cls)

    try:
        table = pd.read_csv(path, sep=separator, compression="infer")
    except pd.errors.EmptyDataError as error:
        raise error_cls(f"{context} file has no parseable content: {path}") from error
    except (pd.errors.ParserError, OSError, UnicodeDecodeError) as error:
        raise error_cls(f"{context} file could not be parsed: {path} ({error})") from error

    if table.shape[1] == 0:
        raise error_cls(f"{context} file contains no columns: {path}")

    if table.shape[0] == 0:
        raise error_cls(f"{context} file contains no data rows: {path}")

    return table


def validate_expression_matrix(
    path: Path,
    gene_id_column: str,
) -> ValidatedExpressionMatrix:
    """Structurally validate an expression matrix.

    Checks performed, per specification Sections 3.1-3.2 and 8.1:
    file existence/readability, presence of the gene-identifier
    column, presence of at least one sample column, no missing gene
    identifiers, no duplicate gene identifiers, no duplicate sample
    column names, and numeric expression values in every sample
    column.

    Args:
        path: Path to the expression matrix file.
        gene_id_column: Name of the column identifying genes.

    Returns:
        A ValidatedExpressionMatrix wrapping the loaded, numerically
        coerced table.

    Raises:
        InvalidExpressionMatrixError: If any structural check fails.
    """
    table = _load_table(path, "Expression matrix", InvalidExpressionMatrixError)

    if gene_id_column not in table.columns:
        raise InvalidExpressionMatrixError(
            f"Expression matrix is missing the configured gene "
            f"identifier column '{gene_id_column}'. Available "
            f"columns: {list(table.columns)}."
        )

    sample_columns = [column for column in table.columns if column != gene_id_column]

    if not sample_columns:
        raise InvalidExpressionMatrixError(
            "Expression matrix contains no sample columns other than "
            f"the gene identifier column '{gene_id_column}'."
        )

    if len(sample_columns) != len(set(sample_columns)):
        seen: set[str] = set()
        duplicates: list[str] = []
        for column in sample_columns:
            if column in seen:
                duplicates.append(column)
            seen.add(column)
        raise InvalidExpressionMatrixError(
            f"Expression matrix contains duplicate sample column "
            f"names: {sorted(set(duplicates))}."
        )

    gene_id_series = table[gene_id_column]

    if gene_id_series.isna().any():
        raise InvalidExpressionMatrixError(
            "Expression matrix contains missing (null) gene "
            f"identifiers in column '{gene_id_column}'."
        )

    duplicated_gene_ids = gene_id_series[gene_id_series.duplicated()].unique().tolist()
    if duplicated_gene_ids:
        raise InvalidExpressionMatrixError(
            f"Expression matrix contains duplicate gene identifiers: "
            f"{sorted(str(gene_id) for gene_id in duplicated_gene_ids)}."
        )

    non_numeric_columns: list[str] = []
    coerced_columns: dict[str, pd.Series] = {}
    for column in sample_columns:
        original = table[column]
        coerced = pd.to_numeric(original, errors="coerce")
        introduced_nan = coerced.isna() & original.notna()
        if introduced_nan.any():
            non_numeric_columns.append(column)
        coerced_columns[column] = coerced

    if non_numeric_columns:
        raise InvalidExpressionMatrixError(
            "Expression matrix contains non-numeric expression values "
            f"in sample column(s): {sorted(non_numeric_columns)}."
        )

    numeric_table = table.copy()
    for column, coerced in coerced_columns.items():
        numeric_table[column] = coerced

    return ValidatedExpressionMatrix(
        file_path=path,
        gene_id_column=gene_id_column,
        sample_columns=tuple(sample_columns),
        gene_ids=tuple(str(gene_id) for gene_id in table[gene_id_column]),
        dataframe=numeric_table,
    )


def validate_metadata(
    path: Path,
    sample_id_column: str,
    group_column: str,
) -> ValidatedMetadata:
    """Structurally validate sample metadata.

    Checks performed, per specification Section 3.3 and 8.1: file
    existence/readability, presence of the sample-identifier and
    group columns, no missing/duplicate sample identifiers, and no
    missing/empty group labels.

    Args:
        path: Path to the metadata file.
        sample_id_column: Name of the column identifying samples.
        group_column: Name of the column identifying each sample's
            biological group.

    Returns:
        A ValidatedMetadata wrapping the loaded table.

    Raises:
        InvalidMetadataError: If any structural check fails.
    """
    table = _load_table(path, "Metadata", InvalidMetadataError)

    for column_name in (sample_id_column, group_column):
        if column_name not in table.columns:
            raise InvalidMetadataError(
                f"Metadata is missing the configured column "
                f"'{column_name}'. Available columns: "
                f"{list(table.columns)}."
            )

    sample_id_series = table[sample_id_column]

    if sample_id_series.isna().any():
        raise InvalidMetadataError(
            "Metadata contains missing (null) sample identifiers in "
            f"column '{sample_id_column}'."
        )

    duplicated_sample_ids = (
        sample_id_series[sample_id_series.duplicated()].unique().tolist()
    )
    if duplicated_sample_ids:
        raise InvalidMetadataError(
            "Metadata contains duplicate sample identifiers: "
            f"{sorted(str(sample_id) for sample_id in duplicated_sample_ids)}."
        )

    group_series = table[group_column]
    empty_group_mask = group_series.isna() | (
        group_series.astype(str).str.strip() == ""
    )

    if empty_group_mask.any():
        affected_samples = (
            table.loc[empty_group_mask, sample_id_column].astype(str).tolist()
        )
        raise InvalidMetadataError(
            "Metadata contains missing or empty group labels for "
            f"sample(s): {sorted(affected_samples)}."
        )

    return ValidatedMetadata(
        file_path=path,
        sample_id_column=sample_id_column,
        group_column=group_column,
        dataframe=table,
    )


def match_samples(
    expression: ValidatedExpressionMatrix,
    metadata: ValidatedMetadata,
) -> SampleMatchResult:
    """Reconcile expression-matrix sample columns against metadata
    sample identifiers.

    Per specification Section 3.4, samples that cannot be reliably
    matched must not be implicitly assigned to a comparison group;
    this function reports (rather than discards) expression-only and
    metadata-only samples, and fails explicitly if there is zero
    usable overlap.

    Args:
        expression: A validated expression matrix.
        metadata: Validated sample metadata.

    Returns:
        A SampleMatchResult with sorted, deterministic tuples of
        matched, expression-only, and metadata-only sample
        identifiers.

    Raises:
        SampleMismatchError: If no sample identifiers overlap between
            the expression matrix and the metadata.
    """
    expression_samples = set(expression.sample_columns)
    metadata_samples = set(
        str(sample_id)
        for sample_id in metadata.dataframe[metadata.sample_id_column]
    )

    matched = sorted(expression_samples & metadata_samples)
    expression_only = sorted(expression_samples - metadata_samples)
    metadata_only = sorted(metadata_samples - expression_samples)

    if not matched:
        raise SampleMismatchError(
            "No samples could be matched between the expression "
            "matrix and the metadata; zero usable sample overlap. "
            f"Expression-only samples: {expression_only}. "
            f"Metadata-only samples: {metadata_only}."
        )

    return SampleMatchResult(
        matched_samples=tuple(matched),
        expression_only_samples=tuple(expression_only),
        metadata_only_samples=tuple(metadata_only),
    )